#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MeshAttributes.h"

void extractPointstoPlot(int NumElx)
{

	int ElNodes[NumElx][2];
	// read in the co-ordinates
	FILE *Nodes = fopen("NodesonXAxis.in","r");
	for (int j = 0; j < NumElx; ++j)
	{
		char nodelist[500];
		int node;
		fgets(nodelist, sizeof(nodelist), Nodes);
		sscanf(nodelist, " %d", &node);
		ElNodes[j][0] = node-1;
	}
	char last_node[500];
	int lastnode;
	fgets(last_node, sizeof(last_node), Nodes);
	sscanf(last_node, "%d", &lastnode);

	fclose(Nodes);

	for (int j = 0; j < NumElx-1; ++j)
	{
		ElNodes[j][1] = ElNodes[j+1][0];
	}	
	ElNodes[NumElx-1][1] = lastnode-1;
	
	int Edges[NumElx];
	for (int i = 0; i < NumElx; ++i)
	{
		int n0 = ElNodes[i][0];
		int n1 = ElNodes[i][1];
		for (int edg = 0; edg < NumEdges; ++edg)
		{
			int e0 = EdgtoVert[edg*2+0];
			int e1 = EdgtoVert[edg*2+1];
			if (((e0==n0) && (e1==n1)) || ((e0==n1) && (e1==n0)))
			{
				Edges[i] = edg;
				break;
			}
	
		}

	}

	int ElEdges[NumElx][2];
	for (int i = 0; i < NumElx; ++i)
	{
		for (int el = 0; el < NumEl; ++el)
		{
		/*	int e0 = EltoEdges[el*3+0];
			int e1 = EltoEdges[el*3+1];
			int e2 = EltoEdges[el*3+2];
*/
			int v0 = EltoVert[el*3+0];
			int v1 = EltoVert[el*3+1];
			int v2 = EltoVert[el*3+2];

			int Ev0 = EdgtoVert[Edges[i]*2];
			int Ev1 = EdgtoVert[Edges[i]*2+1];

			//if (e0 == Edges[i] )
			if ((Ev0 == v0 && Ev1 == v1) || (Ev1==v0 && Ev0 == v1))
			{
				ElEdges[i][0] = el;
				ElEdges[i][1] = 0;
				break;
			}

			//else if (e1 == Edges[i])
			else if ((Ev0 == v1 && Ev1 == v2) || (Ev1 == v1 && Ev0 == v2))
			{
				ElEdges[i][0] = el;
				ElEdges[i][1] = 1;
				break;
			}
			//else if (e2 == Edges[i])
			else if ((Ev0 == v2 && Ev1 == v0) || (Ev1 == v2 && Ev0 == v0))
			{
				ElEdges[i][0] = el;
				ElEdges[i][1] = 2;
				break;
			}
		}
	}

	char filename[] = {"XAxisEdges.out"};
	//sprintf(filename, "extractedEdges_%dElx.txt", NumElx);
	FILE *writefile = fopen(filename,"w");
	for (int i =0; i < NumElx; ++i)
	{
		fprintf(writefile, "%d \t %d \n", ElEdges[i][0], ElEdges[i][1]);

	}
	fclose(writefile);
/*
	FILE *EltoEdgesfile = fopen("EltoEdges.txt","w");
	for (int i =0; i <NumEl; ++i)
	{
		fprintf(EltoEdgesfile, "%d\t%d\t%d\n", EltoEdges[i*3+0], EltoEdges[i*3+1],EltoEdges[i*3+2]);
	}
	fclose(EltoEdgesfile);

	FILE *EdgtoVertfile = fopen("EdgtoVert.txt", "w");
	for (int i = 0; i <NumEdges; ++i)
	{
		fprintf(EdgtoVertfile, "%d\t%d\n", EdgtoVert[i*2+0], EdgtoVert[i*2+1]);
	}
	fclose(EdgtoVertfile);

	FILE *EdgtoElsfile = fopen("EdgtoEls.txt","w");
	for (int i = 0; i < NumEdges; ++i)
	{
		fprintf(EdgtoElsfile, "%d\t%d\n", EdgtoEls[i*2+0], EdgtoEls[i*2+1]);
	}
	fclose(EdgtoElsfile); 
*/
}

