#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <omp.h>
#include <math.h>
#include <mpi.h>
#define NDELAY 600 

#define INFINITE 0x7ff0000000000000

int total_nodes;               	  	//total no. of nodes
int total_iNodes;               	//total no. of input nodes
int total_oNodes;                	//total no. of output nodes
int total_edges;                 	//total no. of edges

typedef struct {
        int node_id;
        int count;
        double iat[20];
} NODE_IAT;

/* Structure of Edge in a graph  */

typedef struct {
        int index_sourceN;         	//source node index in NodeList
        int index_sinkN;            	//sink node index in NodeList
		
        double dei;                     //intrinsic delay of edge
        double den;                     //characteristic delay of edge
        double de;                      //processing delay at edge
}edge;

/* Structure of Node in a graph */

typedef struct node_{
        int index;                      //position of the node is node array
        int fanout;                     //no of edges departing from the node
        int fanin;                      //no of edges sinking to the node

        double dno;                     //intrinsic delay of node
        double dnc;                     //characteristic delay of node
        double m;                       //load
        double dn;                      //processing delay at node


        int  fanout_indexList[20];      //List of fanout edges
        int fanin_indexList[20];        //List of fanin edges
}node;


void find_total_edges_nodes_iNodes_oNodes(char *);
void read_graph_from_input_file(node*  nodeList[],int iNodeList[],int oNodeList[],edge* edgeList[], double at[], double rat[],char *);
double atof(char *);
void calculate_node_delay_of_nodes(node* nodeList[]);
void calculate_edge_delay_of_edges(edge* edgeList[],node * nodeList[]);
void calculate_arrival_time_at_nodes(edge* edgeList[],node * nodeList[],int iNodeList[], double at[]);
void calculate_required_arrival_time_at_nodes(edge* edgeList[],node * nodeList[],int oNodeList[], double rat[]);
void calculate_network_slack_at_nodes(double at[], double rat[], double slack[]);
double delay(double iat);
void master(char *argv[]);
void slave(void);
MPI_Datatype Init_Type_node_MPI (void);
void calculate_no_of_elements_per_process(int n,int p,int dist[],int displs[]);
//void calculate_displacement(int n,int p,int displs[]);

/* Main function */


void main(int argc, char *argv[]) 
{

    // MPI initialization 
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       

    if(!rank) 
      master(argv);
    else 
      slave();

  MPI_Finalize();

}

MPI_Datatype Init_Type_node_MPI () {
     MPI_Datatype NODE_MPI;
     MPI_Datatype type[3] = { MPI_INT, MPI_INT, MPI_DOUBLE };
     int blocklen[3] = { 1, 1, 20 };
     MPI_Aint disp[3];
     disp[0] = 0;
     disp[1] = sizeof(int);
     disp[2] = 2*(sizeof(int));
     MPI_Type_struct (3, blocklen, disp, type, &NODE_MPI);
     return NODE_MPI;
}

void master(char *argv[]) {

        // this work (read the graph,init graph is only for master thread 

        NODE_IAT nodeInfo;   // this is c structure will contain the node package to be sent to slave
        MPI_Datatype NODE_MPI = Init_Type_node_MPI ();
        MPI_Type_commit (&NODE_MPI); 
              
	/*Get the total no. of edges, nodes, iNodes & oNodes from input file */
	find_total_edges_nodes_iNodes_oNodes(argv[1]);
    
	/*Initialize the lists for node, edges, iNode, oNode, at, rat*/

	node** nodeList= (node **) malloc(sizeof(node*)*total_nodes);
	edge** edgeList = (edge **) malloc(sizeof(edge*)*total_edges);;
	int iNodeList[total_iNodes];
	int oNodeList[total_oNodes];
	double *at = (double *) malloc(sizeof(double)*total_nodes);
	double *rat =(double *) malloc(sizeof(double)*total_nodes);
        double start, end;
        double diff_time;

	/*::::::::::::::::INITILAIZE GRAPH FROM INPUT FILE:::::::::::::::::::::::: */

	read_graph_from_input_file(nodeList,iNodeList,oNodeList,edgeList,at,rat,argv[1]);
        printf("Initialized the graph from input file\n"); fflush(stdout);

	start = omp_get_wtime();
	/*::::::::::::::::::::::::: NODE DELAY: dn :::::::::::::::::::::::::::::::: */    
	calculate_node_delay_of_nodes(nodeList);
	printf("Calculated node delay\n"); fflush(stdout);

	/*::::::::::::::Communication/Procesing Delay at Edges:de ::::::::::::::::::*/
	calculate_edge_delay_of_edges(edgeList,nodeList);
	printf("Calculated edge delay\n"); fflush(stdout);
     

        
	/*:::::::::::::::Arrival Time at Nodes - At ::::::::::::::::::::::::::::::::*/
	calculate_arrival_time_at_nodes(edgeList,nodeList,iNodeList,at);
	printf("Calculated arrival time\n"); fflush(stdout);

	/*:::::::::::::: Required Arrival Time at nodes - rat ::::::::::::::::::::::*/
	calculate_required_arrival_time_at_nodes(edgeList,nodeList,oNodeList,rat);
	printf("Calculated required arrival time\n"); fflush(stdout);
     
	/*:::::::::::::::::::::::: Slack at nodes :::::::::::::::::::::::::::::::::::*/

	end = omp_get_wtime();

        diff_time = end-start;
        printf("------------------------------------------------\n");
	printf("It took %lf secs for MPI version.\n", diff_time);
        printf("------------------------------------------------\n");



	//Take argv[1] and remove .in from filename to produce output with .out

	char *output;
	output = strtok( argv[1], "." );
	strcat(output,".out");
	printf( "Output file name is %s \n", output);

	// write outputs to file name with .out extension
	FILE *fp;
	fp=fopen(output, "w");
	int i=0;
	for(i=0;i<total_nodes;i++)
	{
		fprintf(fp,"n%d %.2f %.2f %.2f\n",i,rat[i],at[i],(rat[i]-at[i]));
	}
	fclose(fp);
	printf("Printed the slack in the output file\n"); fflush(stdout);    
}


void calculate_network_slack_at_nodes(double at[], double rat[], double slack[])
{
	int i = 0;
	for (i =0; i < total_nodes; i++)
	{
		slack[i]= rat[i] - at[i];
	}
	
}


void calculate_arrival_time_at_nodes(edge* edgeList[],node* nodeList[], int iNodeList[], double at[])
{
    node** nodeQueue = (node **)malloc(sizeof(node*)*total_nodes);
    int front_indexQ_temp,front_indexQ=0;
    int rear_indexQ_temp,rear_indexQ=0; 
    int *fanin_nodes = (int *)malloc(sizeof(int)*total_nodes);     // fanin for any node 
    int i=0,j=0;
    int index_sinkN=0, index_sourceN=0;	
    double at_temp;

    /*--------------------MPI ---------------------*/
    int nprocs,rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       
    MPI_Datatype NODE_MPI = Init_Type_node_MPI ();
    MPI_Type_commit (&NODE_MPI);
    int element_dist[nprocs];
    int displs[nprocs];
    /*----------------------------------------------*/


    for(i=0;i<total_nodes;i++)
    {
       fanin_nodes[i]= nodeList[i]->fanin;
    }

    // insert in-nodes in Queue 
    for (i=0;i<total_iNodes;i++)
    {
       nodeQueue[rear_indexQ] = nodeList[iNodeList[i]]; 
       rear_indexQ++;
    }

    while(rear_indexQ > front_indexQ)
    {

          front_indexQ_temp = front_indexQ;
          rear_indexQ_temp = rear_indexQ;
           
          int elements_on_level = (rear_indexQ - front_indexQ);

	  // add who will get whom and displacement 
          calculate_no_of_elements_per_process(elements_on_level,nprocs,element_dist,displs);             
		
          NODE_IAT* node_send_PACKAGE = (NODE_IAT*) malloc (sizeof(NODE_IAT)*elements_on_level);
          NODE_IAT* node_recv_PACKAGE = (NODE_IAT*) malloc (sizeof(NODE_IAT)*(element_dist[rank]));   // to be changed
          double* node_send_PACKAGE_G = (double *) malloc (sizeof(double)*(element_dist[rank]));
          double* node_recv_PACKAGE_G = (double *) malloc (sizeof(double)*elements_on_level);
 
          int mn=0;
          for (i=front_indexQ_temp;i<rear_indexQ_temp;i++) {           
            int not_on_level0=0;
            node_send_PACKAGE[mn].node_id = nodeQueue[i]->index;    // Filling the node index in package
            node_send_PACKAGE[mn].count = nodeQueue[i]->fanin;         // Filling the no. of fanin in package
            
            for(j=0;j<nodeQueue[i]->fanin;j++) {
                at_temp = 0;
                index_sourceN = edgeList[((nodeQueue[i])->fanin_indexList[j])]->index_sourceN;
                index_sinkN = nodeQueue[i]->index;
                at_temp = at[index_sourceN] + nodeList[index_sourceN]->dn+
                             edgeList[((nodeQueue[i])->fanin_indexList[j])]->de;
                node_send_PACKAGE[mn].iat[j] = at_temp;              // Filling the values from node parents   
                not_on_level0 =1;
            }
            if(!not_on_level0) { node_send_PACKAGE[mn].iat[0] = at[nodeQueue[i]->index]; node_send_PACKAGE[mn].count =1; }
            mn++;
          }
         //broadcasting the total number of element on a level
         MPI_Bcast(&elements_on_level, 1, MPI_INT, 0, MPI_COMM_WORLD);

         //scattering the elements on a level 
         MPI_Scatterv(node_send_PACKAGE,element_dist,displs, NODE_MPI, node_recv_PACKAGE,element_dist[rank], NODE_MPI, 0, MPI_COMM_WORLD);

         for (i=0;i<element_dist[rank];i++) {
           int not_on_level0=0;
           double at_max =0;
           for (j=0;j<node_recv_PACKAGE[i].count;j++) {
               at_max = (at_max > (node_recv_PACKAGE[i].iat[j]))?at_max:(node_recv_PACKAGE[i].iat[j]);
           }
           at_max = delay(at_max);
           node_send_PACKAGE_G[i]=at_max;
         }

         //gather here 
         MPI_Gatherv(node_send_PACKAGE_G,element_dist[rank],MPI_DOUBLE, node_recv_PACKAGE_G,element_dist,displs,MPI_DOUBLE, 0, MPI_COMM_WORLD);
         mn=0; 
         for (i=front_indexQ_temp;i<rear_indexQ_temp;i++) {
             at[nodeQueue[i]->index] = node_recv_PACKAGE_G[mn++];
         }

         for (i=front_indexQ_temp;i<rear_indexQ_temp;i++) {

            for(j=0;j<nodeQueue[i]->fanout;j++) {
                index_sinkN =edgeList[((nodeQueue[i])->fanout_indexList[j])]->index_sinkN;
                if(--fanin_nodes[index_sinkN]==0)
                   nodeQueue[rear_indexQ++]=nodeList[index_sinkN];
            }
          }
          front_indexQ = rear_indexQ_temp;

         free(node_send_PACKAGE);
         free(node_recv_PACKAGE);
         free(node_send_PACKAGE_G);
         free(node_recv_PACKAGE_G);
        
    }
        free(fanin_nodes);
        //broadcasting -1 to end the loop
        int quit = -1;
        MPI_Bcast(&quit, 1, MPI_INT, 0, MPI_COMM_WORLD);
}



void calculate_required_arrival_time_at_nodes(edge* edgeList[],node * nodeList[],int oNodeList[], double rat[])
{
    node** nodeQueue = (node **)malloc(sizeof(node *)*total_nodes);
    int front_indexQ_temp,front_indexQ=0;
    int rear_indexQ_temp, rear_indexQ=0;
    int *fanout_nodes =(int *)malloc(sizeof(int)*total_nodes);     // fanout for any node
    int i=0,j=0;
    int index_sinkN=0, index_sourceN=0;
    double rat_temp;

    /*--------------------MPI ---------------------*/
    int nprocs,rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype NODE_MPI = Init_Type_node_MPI ();
    MPI_Type_commit (&NODE_MPI);
    int element_dist[nprocs];
    int displs[nprocs];
    /*----------------------------------------------*/

    for(i=0;i<total_nodes;i++)
    {
       fanout_nodes[i]= nodeList[i]->fanout;
    }

    // insert in-nodes in Queue
    for (i=0;i<total_oNodes;i++)
    {
        nodeQueue[rear_indexQ] = nodeList[oNodeList[i]];
        rear_indexQ++;
    }

    while(rear_indexQ > front_indexQ)
    {
          front_indexQ_temp = front_indexQ;
          rear_indexQ_temp = rear_indexQ;         

          int elements_on_level = (rear_indexQ - front_indexQ);

          // add who will get whom and displacement 
          calculate_no_of_elements_per_process(elements_on_level,nprocs,element_dist,displs);

          NODE_IAT* node_send_PACKAGE = (NODE_IAT*) malloc (sizeof(NODE_IAT)*elements_on_level);
          NODE_IAT* node_recv_PACKAGE = (NODE_IAT*) malloc (sizeof(NODE_IAT)*(element_dist[rank]));   // to be changed
          double* node_send_PACKAGE_G = (double *) malloc (sizeof(double)*(element_dist[rank]));
          double* node_recv_PACKAGE_G = (double *) malloc (sizeof(double)*elements_on_level);


          int mn=0;
          for (i=front_indexQ_temp;i<rear_indexQ_temp;i++) {
            int not_on_level0=0;
            node_send_PACKAGE[mn].node_id = nodeQueue[i]->index;    // Filling the node index in package
            node_send_PACKAGE[mn].count = nodeQueue[i]->fanout;         // Filling the no. of fanout in package

            for(j=0;j<nodeQueue[i]->fanout;j++) {
                rat_temp = 0;
                index_sourceN = edgeList[((nodeQueue[i])->fanout_indexList[j])]->index_sinkN;
                index_sinkN = nodeQueue[i]->index;
                rat_temp = (rat[index_sourceN] -
                                       (nodeList[index_sinkN]->dn +
                                       edgeList[((nodeQueue[i])->fanout_indexList[j])]->de));
                node_send_PACKAGE[mn].iat[j] = rat_temp;              // Filling the values from node parents   
                not_on_level0 =1;
            }
            if(!not_on_level0) { node_send_PACKAGE[mn].iat[0] = rat[nodeQueue[i]->index]; node_send_PACKAGE[mn].count =1; }
            mn++;
          }


         //broadcasting the total number of element on a level
         MPI_Bcast(&elements_on_level, 1, MPI_INT, 0, MPI_COMM_WORLD);

         //scattering the elements on a level 
         MPI_Scatterv(node_send_PACKAGE,element_dist,displs, NODE_MPI, node_recv_PACKAGE,element_dist[rank], NODE_MPI, 0, MPI_COMM_WORLD);

         for (i=0;i<element_dist[rank];i++) {
           int not_on_level0=0;
           double rat_min =INFINITE;
           for (j=0;j<node_recv_PACKAGE[i].count;j++) {
               rat_min = (rat_min < (node_recv_PACKAGE[i].iat[j]))?rat_min:(node_recv_PACKAGE[i].iat[j]);
           }
           //rat_min = delay(rat_min);
           node_send_PACKAGE_G[i]=rat_min;
         }

         //gather here 
         MPI_Gatherv(node_send_PACKAGE_G,element_dist[rank],MPI_DOUBLE, node_recv_PACKAGE_G,element_dist,displs,MPI_DOUBLE, 0, MPI_COMM_WORLD);
         mn=0;
         for (i=front_indexQ_temp;i<rear_indexQ_temp;i++) {
             rat[nodeQueue[i]->index] = node_recv_PACKAGE_G[mn++];
         }

         for (i= front_indexQ_temp;i<rear_indexQ_temp;i++) {

            for(j=0;j<nodeQueue[i]->fanin;j++) {
                index_sinkN =edgeList[((nodeQueue[i])->fanin_indexList[j])]->index_sourceN;
                if(--fanout_nodes[index_sinkN]==0)
                   nodeQueue[rear_indexQ++]=nodeList[index_sinkN];
            }
        }
        front_indexQ = rear_indexQ_temp;

        free(node_send_PACKAGE);
        free(node_recv_PACKAGE);
        free(node_send_PACKAGE_G);
        free(node_recv_PACKAGE_G);
   }
        free(fanout_nodes);
        // Send Bcast to slaves to Quit 
        int quit = -2;
        MPI_Bcast(&quit, 1, MPI_INT, 0, MPI_COMM_WORLD);

}



void calculate_node_delay_of_nodes(node*  nodeList[])
{
	int i;
	for (i=0;i<total_nodes;i++)
	{
		nodeList[i]->dn = nodeList[i]->dno + 
				  (nodeList[i]->fanout)*(nodeList[i]->dnc);	
		
        }


}

/*
Input: List of Nodes, List of Edges
Output: Node delay at Edges
Formula: de = dei + m den
*/



void calculate_edge_delay_of_edges(edge* edgeList[],node * nodeList[])
{
	int i;
	double m; // Load of Sink Noad
	for (i=0;i<total_edges;i++)
	{
		m = nodeList[(edgeList[i]->index_sinkN)]->m;
		edgeList[i]->de = edgeList[i]->dei + (m*(edgeList[i]->den));

	}

}



void read_graph_from_input_file(node*  nodeList[],int iNodeList[],int oNodeList[],edge* edgeList[],double at[], double rat[], char *fName)
{
  FILE *f = fopen(fName, "r" );
  char line[128];
  int iNode_index =0;
  int oNode_index =0;
  int node_index =0;
  int edge_index =0;
  int fanout =0;
  int fanin=0;
  double temp;
  int i1,i2,i3;
  char node_name[20];
  char node_name2[20];
  double dno,dnc,dei,den,at_tmp,rat_tmp,m; 
  char type;
  
  for(i3=0;i3<total_nodes;i3++)
  {
      rat[i3]=INFINITE;
  }

  while(fgets(line,128,f)!=NULL)
  {
     switch(line[0])
     {
        case 'i':
                sscanf(line,"%c %s",&type,node_name);
		iNodeList[iNode_index] = atoi(node_name+1);
		iNode_index++;
                break;
        case 'o':
                sscanf(line,"%c %s",&type,node_name);
		oNodeList[oNode_index] = atoi(node_name+1);
                oNode_index++;		
                break;
        case 'n':
                sscanf(line,"%c %s %lf %lf %lf",&type,node_name,&dno,&dnc,&m);
                int indx = atoi(node_name+1);
                nodeList[indx]=(node *)malloc(sizeof(node));
		(nodeList[indx])->index = indx;
			
		nodeList[indx]->dno = dno;

		nodeList[indx]->dnc = dnc;

		nodeList[indx]->m = m;
                
                nodeList[indx]->fanout =0;

                nodeList[indx]->fanin = 0;
                break;
        case 'e':
		sscanf(line,"%c %s %s %lf %lf",&type,node_name,node_name2,&dei,&den);
		edgeList[edge_index]=(edge *)malloc(sizeof(edge));		

		i1 = atoi(node_name+1);
		edgeList[edge_index]->index_sourceN = i1;
		fanout = nodeList[i1]->fanout;
		nodeList[i1]->fanout_indexList[fanout]=edge_index;
		nodeList[i1]->fanout++;

		i1 = atoi(node_name2+1);
                edgeList[edge_index]->index_sinkN = i1;
                fanin = nodeList[i1]->fanin;
                nodeList[i1]->fanin_indexList[fanin]=edge_index;
                nodeList[i1]->fanin++;

		edgeList[edge_index]->dei = dei;

		edgeList[edge_index]->den = den;

		edge_index++;

                break;
        case 'a':
                sscanf(line,"%c %s %lf",&type,node_name,&at_tmp);
                i2 = atoi(node_name+1);
		at[i2] = at_tmp;
                break;
        case 'r':
		sscanf(line,"%c %s %lf",&type,node_name,&rat_tmp);
                i2 = atoi(node_name+1);	
		rat[i2]=rat_tmp;
                break;
        default:
                if(!total_nodes) { 
               /*
                total_nodes = atoi(buf);
                fscanf(f, "%s", buf);
                total_iNodes = atoi(buf);
                fscanf(f, "%s", buf);
                total_oNodes = atoi(buf);
		*/
                }
                break;
     }
  }
  fclose(f);
}

void find_total_edges_nodes_iNodes_oNodes(char *fName)
{
	FILE *file = fopen(fName, "r");
	char line [ 128 ]; /* or other suitable maximum line size */
	int i=0;
	while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */
	{		
		if(line[0]=='n') {total_nodes++;}
		else if(line[0]=='e') {total_edges++;}
		else if(line[0]=='i') {total_iNodes++;}
		else if(line[0]=='o') {total_oNodes++;}
	}

	fclose ( file );

}


double delay(double iat) {
  int i;
  double sc = 0.0;

  for(i = 0; i < NDELAY; i++) {
    sc += 1.0 / (pow(tan(i * M_PI/NDELAY), 2.0) + 1.0) + pow(sin(i * M_PI/NDELAY), 2.0);
  }
  sc /= NDELAY;

  return(iat*sc);
}


void calculate_no_of_elements_per_process(int n,int p,int dist[],int displs[]) {
  int i;
  int r = n%p;
  
  for (i= 0;i<p;i++) {
     if(r==0) dist[i] = (int)(n/p);
     else if(i<r) dist[i] = (int)(n/p) +1;
     else dist[i] = (int)(n/p);
  }
  displs[0]=0;
  int sum=0;
  for (i=0;i<(p-1);i++) {
     sum = sum + dist[i];
     displs[i+1]=sum;
  }
}



void slave() {
   
        int nprocs,rank;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Datatype NODE_MPI = Init_Type_node_MPI ();
        MPI_Type_commit (&NODE_MPI);
        int element_dist[nprocs];
        int displs[nprocs];
        int elements_on_level;
        int i,j;
        int rat_turn=0;

   while(1) {
        // recieve broadcast from master 
        MPI_Bcast(&elements_on_level, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(elements_on_level == -1) { 
	  rat_turn=1;  // it's rat's turn now 
          MPI_Bcast(&elements_on_level, 1, MPI_INT, 0, MPI_COMM_WORLD);
        } 
         
        if(elements_on_level == -2) { break; } // rat is done, quit 

        calculate_no_of_elements_per_process(elements_on_level,nprocs,element_dist,displs);

        NODE_IAT* node_send_PACKAGE = NULL;
        NODE_IAT* node_recv_PACKAGE = (NODE_IAT*) malloc (sizeof(NODE_IAT)*(element_dist[rank]));  
        double* node_send_PACKAGE_G = (double *) malloc (sizeof(double)*(element_dist[rank]));
        double* node_recv_PACKAGE_G = (double *) malloc (sizeof(double)*elements_on_level);

        
        MPI_Scatterv(node_send_PACKAGE,element_dist,displs, NODE_MPI, node_recv_PACKAGE,element_dist[rank], NODE_MPI, 0, MPI_COMM_WORLD);
    
        if(rat_turn) {
           for (i=0;i<element_dist[rank];i++) {
               double rat_min =INFINITE;
               for (j=0;j<node_recv_PACKAGE[i].count;j++) {
                   rat_min = (rat_min < (node_recv_PACKAGE[i].iat[j]))?rat_min:(node_recv_PACKAGE[i].iat[j]);
               }
              rat_min = delay(rat_min);
              node_send_PACKAGE_G[i]=rat_min;
           }
        }
        else {
           for (i=0;i<element_dist[rank];i++) {
               double at_max =0;
               for (j=0;j<node_recv_PACKAGE[i].count;j++) {
                   at_max = (at_max > (node_recv_PACKAGE[i].iat[j]))?at_max:(node_recv_PACKAGE[i].iat[j]);
               }
               at_max = delay(at_max);
               node_send_PACKAGE_G[i]=at_max;
           }
        }
       MPI_Gatherv(node_send_PACKAGE_G,element_dist[rank],MPI_DOUBLE, node_recv_PACKAGE_G,element_dist,displs,MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
}
