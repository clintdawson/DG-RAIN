#include <stdio.h>
#include <stdlib.h>
#include <math.h>	
#include "mathfunctions.h"
#include "MeshAttributes.h"

int NumEl;
int NumNodes;
int NumEdges;
double *x_coord;
double *y_coord;
double *z;
double *nFriction;
int *EltoVert;
int *EdgtoEls;
int *EdgtoVert;
int *BdryPrescribed;
int numInflowEdges;
int numOutflowEdges;
int secondSetInflowEdges;
int *InflowBdryNodes1;
int *InflowBdryNodes2;
int *OutflowBdryNodes;

double minEdgLength = 9999999;
double max_lambda;

void *xcalloc(int items, int size)
{
	void *ptr = calloc(items, size);
	if (ptr == NULL)
	{
		printf("Unable to allocate memory\n");
		exit(EXIT_FAILURE);
	}
	return ptr;
}

void store_mesh(char *FullMesh)
{
	// read in the co-ordinates
	FILE *Coordinates = fopen(FullMesh,"r");
	
	// throw away the first line of the file
	char firstline[1000];
	fgets(firstline, sizeof(firstline), Coordinates);

	char secondline[1000];
	fgets(secondline, sizeof(secondline), Coordinates);
	sscanf(secondline, "%d  %d", &NumEl, &NumNodes);	
	
	// array to store the x, y and z values at the nodes
	x_coord = xcalloc(NumNodes,sizeof(double));
	y_coord = xcalloc(NumNodes,sizeof(double));
	z = xcalloc(NumNodes,sizeof(double));
	nFriction = xcalloc(NumNodes,sizeof(double));

	// variable to temporarily store the node numbers while reading the file
	int nodes;
	
	int j;
	for (j=0; j<NumNodes; ++j)
	{	
		char grid[1000];
		fgets(grid, sizeof(grid), Coordinates);
		sscanf(grid, "  %d %lf %lf %lf %lf", &nodes, &x_coord[j], &y_coord[j], &z[j], &nFriction[j]);
		z[j] = z[j];

	}

	EltoVert = xcalloc(3*NumEl,sizeof(int));
	int numNodesPerEl;
	int ElNum;
	for (j=0; j < NumEl; ++j)
	{
		char conn_table[1000];
		fgets(conn_table, sizeof(conn_table), Coordinates);
		int v1, v2, v3;
		sscanf(conn_table, "%d %d %d %d %d", &ElNum, &numNodesPerEl, &v1, &v2, &v3);
		// switch to indexing from 0
		EltoVert[j*3+0] = v1-1;
		EltoVert[j*3+1] = v2-1;
		EltoVert[j*3+2] = v3-1;
	}

	// throw away the next line of the file
	char throwawayline[1000];
	fgets(throwawayline, sizeof(throwawayline), Coordinates);

	// get the number of first inflow boundary edges
	char newline[1000];
	fgets(newline, sizeof(newline), Coordinates);
	sscanf(newline, "%d", &numInflowEdges);
	InflowBdryNodes1 = xcalloc(2*numInflowEdges,sizeof(int));

	for (int k = 0; k < numInflowEdges; ++k)
	{
		char nodenums[1000];
		fgets(nodenums, sizeof(nodenums), Coordinates);
		int num1, num2;
		sscanf(nodenums, "%d %d", &num1, &num2);
		InflowBdryNodes1[2*k] = num1-1;
		InflowBdryNodes1[2*k+1] = num2-1;
	}	

	// get the number of second set of inflow edges
	char newinflowline[1000];
	fgets(newinflowline, sizeof(newinflowline), Coordinates);
	sscanf(newinflowline, "%d", &secondSetInflowEdges);
	InflowBdryNodes2 = xcalloc(2*secondSetInflowEdges,sizeof(int));

	for (int k = 0; k < secondSetInflowEdges; ++k)
	{
		char nodenums[1000];
		fgets(nodenums, sizeof(nodenums), Coordinates);
		int num1, num2;
		sscanf(nodenums, "%d %d", &num1, &num2);
		InflowBdryNodes2[2*k] = num1-1;
		InflowBdryNodes2[2*k+1] = num2-1;
	}	
	
	// get the number of outflow edges
	char outflowline[1000];
	fgets(outflowline, sizeof(outflowline), Coordinates);
	sscanf(outflowline, "%d", &numOutflowEdges);
	OutflowBdryNodes = xcalloc(2*numOutflowEdges,sizeof(int));

	for (int k = 0; k < numOutflowEdges; ++k)
	{
		char nodenums[1000];
		fgets(nodenums, sizeof(nodenums), Coordinates);
		int num1, num2;
		sscanf(nodenums, "%d %d", &num1, &num2);
		OutflowBdryNodes[2*k] = num1-1;
		OutflowBdryNodes[2*k+1] = num2-1;
	}


	fclose(Coordinates);

}

void extract_mesh_attributes()
{
	// For each node, count the number of elements that share this node

//	int NtoNodes[NumNodes][10];
//	int NodesCount[NumNodes];

		int **NtoNodes = xcalloc(NumNodes,sizeof(int*));
		for (int s = 0; s < NumNodes; s++)
		{
			NtoNodes[s] = xcalloc(10,sizeof(int));
		}
		int *NodesCount = xcalloc(NumNodes,sizeof(int));
	
{	
	//int ElCount[NumNodes];	
	int *ElCount = xcalloc(NumNodes,sizeof(int));
	for (int n = 0; n < NumNodes; ++n)
	{
		//ElCount[n] = 0;
		for (int j = 0; j < 10; ++j)
			NtoNodes[n][j] = -1;
	}
	for (int el = 0; el<NumEl; ++el)
	{
		for (int n = 0; n < 3; ++n)
		{
			int node = EltoVert[el*3+n];
			ElCount[node] = ElCount[node]+1;
			int numConnNodes = 0;
			while(NtoNodes[node][numConnNodes] >= 0)
			{
				++numConnNodes;
			}
			
			int n1 = (n+1)%3;
			int n2 = (n+2)%3;
			int neighb0 = EltoVert[el*3+n1];
			int neighb1 = EltoVert[el*3+n2];
			int neighbfound = 0;
			int foundneighb[2] = {-1,-1};
			for (int j = 0; j < 10; ++j)
			{
				if (NtoNodes[node][j] < 0 || neighbfound == 2)
					break;
				if (NtoNodes[node][j] == neighb0)
				{
					neighbfound+=1;
					foundneighb[0] = 1;
				}
				if (NtoNodes[node][j] == neighb1)
				{
					neighbfound+=1;
					foundneighb[1] = 1;
				}
				
			}
			if (foundneighb[0] == -1)
			{
				NtoNodes[node][numConnNodes] = neighb0;
				++numConnNodes;
			}
			if (foundneighb[1] == -1)
			{
				NtoNodes[node][numConnNodes] = neighb1;
				++numConnNodes;
			}
			NodesCount[node] = numConnNodes;

		}

	}

	free(ElCount);

}	// Euler's relation: NumNodes - NumEdges + Numfaces = 1 -> NumEdges = NumNodes + NumEl - 1;
	NumEdges = NumNodes + NumEl - 1;
	EdgtoEls = xcalloc(2*NumEdges,sizeof(int));
	EdgtoVert = xcalloc(2*NumEdges,sizeof(int));
	BdryPrescribed = xcalloc(NumEdges,sizeof(int));
	

	int EdgCount = 0;
	for (int n = 0; n < NumNodes; ++n)
	{
		int NumConnEdges = NodesCount[n];
		//if (NumConnEdges > 6)
		//	printf("NumConnEdges = %d \n", NumConnEdges);
		for (int k = 0; k < NumConnEdges; ++k)
		{
			int connNode = NtoNodes[n][k];
			if (connNode > n)
			{
				EdgtoVert[EdgCount*2] = n;
				EdgtoVert[EdgCount*2+1] = connNode;
				EdgCount += 1;

			}
		//	if (NumConnEdges > 6)
		//	{
		//		printf("nodenum = %d \t conn_node = %d\n", n, connNode);
		//	}
	
		}
	
		//if ( NumConnEdges > 6)
		//	exit(1);
	}

	if (EdgCount != NumEdges)	
	{
		printf("EdgCount = %d\n", EdgCount);
		printf("Total number of edges is not correct\n");
		exit(EXIT_FAILURE);
	}	

	// Create EdgtoEls
	for (int e = 0; e < NumEdges; ++e)
	{
		EdgtoEls[e*2] = -1;
		EdgtoEls[e*2+1] = -1;
	}
	
	for (int e = 0; e < NumEdges; ++e)
	{
		int v0 = EdgtoVert[e*2];
		int v1 = EdgtoVert[e*2+1];
		int elfound = 0;
		for (int el = 0; el < NumEl; ++el)
		{
			for (int k = 0; k < 3; ++k)
			{
				int n0 = EltoVert[el*3+k];
				int ind1 = (k+1)%3;
				int n1 = EltoVert[el*3+ind1];
				if (((n0 == v0) && (n1 == v1))|| ((n0 == v1) &&(n1==v0)))
				{
					if (EdgtoEls[e*2] < 0)
						EdgtoEls[e*2] = el;
					else if (EdgtoEls[e*2+1] < 0)
						EdgtoEls[e*2+1] = el;
					elfound+=1;
					break;
				}
			}
			if (elfound==2)
				break;
			if ((el == (NumEl-1)) && elfound == 1)
				EdgtoEls[e*2+1] = EdgtoEls[e*2];

		 }

	}

	for (int e = 0; e < NumEdges; ++e)
	{
		int n0 = EdgtoVert[e*2+0];
		int n1 = EdgtoVert[e*2+1];

		double edgLength = sqrt((x_coord[n1]-x_coord[n0])*(x_coord[n1]-x_coord[n0])+(y_coord[n1]-y_coord[n0])*(y_coord[n1]-y_coord[n0]));
		minEdgLength = min(minEdgLength,edgLength);

		int InflowBoundary = 0;
		int OutflowBoundary = 0;
		for (int k = 0 ; k < numInflowEdges; ++k)
		{
			int b0 = InflowBdryNodes1[2*k];
			int b1 = InflowBdryNodes1[2*k+1];
			if (((b0==n0) && (b1 ==n1)) || ((b0==n1) && (b1==n0)))
			{
				InflowBoundary = 1;
				break;
			}

		}

		for (int k = 0; k < secondSetInflowEdges; ++k)
		{
			int b0 = InflowBdryNodes2[2*k];
			int b1 = InflowBdryNodes2[2*k+1];
			if (((b0==n0) && (b1 ==n1)) || ((b0==n1) && (b1==n0)))
			{
				InflowBoundary = 2;
				break;
			}
		}

		if (!InflowBoundary)
		{
			for (int k = 0; k < numOutflowEdges; ++k)
			{
				int b0 = OutflowBdryNodes[2*k];
				int b1 = OutflowBdryNodes[2*k+1];
				if (((b0==n0) && (b1 ==n1)) || ((b0==n1) && (b1==n0)))
				{
					OutflowBoundary = 1;
					break;
				}

			}
		}

		if (InflowBoundary == 1)
			BdryPrescribed[e] = 1;
		else if(InflowBoundary == 2)
			BdryPrescribed[e] = 3;
		else if (OutflowBoundary)
			BdryPrescribed[e] = 2;
		else
			BdryPrescribed[e] = 0;

	}

	free(NtoNodes);
	free(NodesCount);
	free(InflowBdryNodes1);
	free(InflowBdryNodes2);
	free(OutflowBdryNodes);
}
