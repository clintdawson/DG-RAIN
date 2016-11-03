#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ChannelsAndJunctions.h"
#include "MeshAttributes.h"

/*********************************************************************************************************//**
* @file create_channel_network.c
* 
* This file contains code to create a channel network from the grid files provided for the channels and the junctions. 
* *********************************************************************************************************/


/*******************************************************************************************************//**
* Function that reads in the grid files and creates a channel structure for each %channel in the network 
* and a junction structure for each %junction in the network. The channel and the junction structures 
* contain all the information about channels and the junctions respectively and their connectivity. 
* The channel structures are then stored in an array called ChannelList and the junction structures are
* stored in an array called JunctionList.
* @param [in] ChannelNodes string name of the gid file for channels
* @param [in] JunctionNodes string name of the grid file for junctions
*
* ********************************************************************************************************/

void create_channel_network(char *ChannelNodes, char *JunctionNodes)
{
	FILE *ChanNodes = fopen(ChannelNodes, "r");
	char line[100];
	fgets(line, sizeof(line),ChanNodes);	
	sscanf(line, "%d", &NumChannels);

	// allocate space for an array to store channels
	ChannelList = malloc(NumChannels*sizeof(struct channel*));

	for (int i = 0; i < NumChannels; ++i)
	{
		char newline[100];
		int NumNodes;
		fgets(newline, sizeof(newline), ChanNodes);
		sscanf(newline, "%d", &NumNodes);

		double x[NumNodes], y[NumNodes], z[NumNodes], b[NumNodes], n[NumNodes];
	
		for (int j = 0; j < NumNodes; ++j)
		{
			char grid[500];
			fgets(grid, sizeof(grid), ChanNodes);
			sscanf(grid, "%lf %lf %lf %lf %lf", &x[j], &y[j], &z[j], &b[j], &n[j]);

		}

		ChannelList[i] = malloc(sizeof(struct channel));

		ChannelList[i]->NumNodes = NumNodes;
		ChannelList[i]->NumEl = NumNodes-1;

		ChannelList[i]->x = malloc(NumNodes*sizeof(double));
		ChannelList[i]->y = malloc(NumNodes*sizeof(double));
		ChannelList[i]->z = malloc(NumNodes*sizeof(double));
		ChannelList[i]->b = malloc(NumNodes*sizeof(double));
		ChannelList[i]->nFriction = malloc(NumNodes*sizeof(double));
		
		for (int j = 0; j < NumNodes; ++j)
		{
			ChannelList[i]->x[j] = x[j];
			ChannelList[i]->y[j] = y[j];
			ChannelList[i]->z[j] = z[j];
			ChannelList[i]->b[j] = b[j];
			ChannelList[i]->nFriction[j] = n[j];
		}

		double mindh = 9999999999;
		for (int j = 0; j < NumNodes-1; ++j)
		{
			double dh = sqrt((x[j+1]-x[j])*(x[j+1]-x[j])+(y[j+1]-y[j])*(y[j+1]-y[j]));
			mindh = fmin(mindh, dh);
		}	
	
		ChannelList[i]->mindh = mindh;
	} // end loop for Channels

	fclose(ChanNodes);

	FILE *JuncNodes = fopen(JunctionNodes, "r");
	char juncLine[100];
	fgets(juncLine, sizeof(line), JuncNodes);
	sscanf(juncLine, "%d", &NumJunctions);

	// allocate space for an array to store junctions
	JunctionList = malloc(NumJunctions*sizeof(struct junction*));
	for (int i = 0; i < NumJunctions; ++i)
	{
		char newline[500];
		int NumEl, NumNodes;
		fgets(newline, sizeof(newline), JuncNodes);
		sscanf(newline, "%d\t%d", &NumEl, &NumNodes);
		
		double x[NumNodes], y[NumNodes], z[NumNodes], n[NumNodes];

		// read and store the coordinates
		for (int j = 0; j < NumNodes; ++j)
		{
			char grid1[500];
			int nodeNum;
			fgets(grid1, sizeof(grid1), JuncNodes);
			sscanf(grid1, "%d %lf %lf %lf %lf", &nodeNum, &x[j], &y[j], &z[j], &n[j]);
		}

		JunctionList[i] = malloc(sizeof(struct junction));

		JunctionList[i]->NumNodes = NumNodes;
		JunctionList[i]->NumEl = NumEl;

		JunctionList[i]->x = malloc(NumNodes*sizeof(double));
		JunctionList[i]->y = malloc(NumNodes*sizeof(double));
		JunctionList[i]->z = malloc(NumNodes*sizeof(double));
		JunctionList[i]->nFriction = malloc(NumNodes*sizeof(double));

		for (int j = 0; j < NumNodes; ++j)
		{
			JunctionList[i]->x[j] = x[j];
			JunctionList[i]->y[j] = y[j];
			JunctionList[i]->z[j] = z[j];
			JunctionList[i]->nFriction[j] = n[j];
		}

		// read and store the adjacency matrix
		JunctionList[i]->EltoVert = malloc(3*NumEl*sizeof(int));
		int EtoV[NumEl][3];
		for (int j = 0; j < NumEl; ++j)
		{
			char table[100];
			int elNum, NumNodesPerEl;
			fgets(table, sizeof(table), JuncNodes);
			sscanf(table, "%d %d %d %d %d", &elNum, &NumNodesPerEl, &EtoV[j][0], &EtoV[j][1], &EtoV[j][2]);
		}
	
		// Euler's relation: NumNodes - NumEdges + Numfaces = 1 -> NumEdges = NumNodes + NumEl - 1;
		int NumEdges = NumNodes + NumEl - 1;
		JunctionList[i]->TotalNumEdges = NumEdges;
		JunctionList[i]->EdgtoEls = malloc(2*NumEdges*sizeof(int));
		JunctionList[i]->EdgtoVert = malloc(2*NumEdges*sizeof(int));

		int NumInflowOutflow;
		char inOut[100];
		fgets(inOut, sizeof(inOut), JuncNodes);
		sscanf(inOut, "%d", &NumInflowOutflow);
		int InOutNodes[NumInflowOutflow][4];

		for (int j = 0; j < NumInflowOutflow; ++j)
		{
			char bdry[500];
			fgets(bdry, sizeof(bdry), JuncNodes);
			sscanf(bdry, "%d %d %d %d", &InOutNodes[j][0], &InOutNodes[j][1], &InOutNodes[j][2], &InOutNodes[j][3]); 
			InOutNodes[j][0] -=1; 
			InOutNodes[j][1] -=1; 
		}

	
		for (int j = 0; j < NumEl; ++j)
		{
			for (int k = 0; k < 3; ++k)		
			{
				JunctionList[i]->EltoVert[j*3+k] = EtoV[j][k]-1;
			}
		}		
		
		// create EdgtoVert and EltoEdges	
		int NtoNodes[NumNodes][20];
		int NodesCount[NumNodes];
	
		int ElCount[NumNodes];	
		for (int n = 0; n < NumNodes; ++n)
		{
			ElCount[n] = 0;
			for (int j = 0; j < 10; ++j)
				NtoNodes[n][j] = -1;
		}
		for (int el = 0; el<NumEl; ++el)
		{
			for (int n = 0; n < 3; ++n)
			{
				int node = JunctionList[i]->EltoVert[el*3+n];
				ElCount[node] = ElCount[node]+1;
				int numConnNodes = 0;
				while(NtoNodes[node][numConnNodes] >= 0)
				{
					++numConnNodes;
				}
				
				int n1 = (n+1)%3;
				int n2 = (n+2)%3;
				int neighb0 = JunctionList[i]->EltoVert[el*3+n1];
				int neighb1 = JunctionList[i]->EltoVert[el*3+n2];
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
		
			
		int EdgCount = 0;
		for (int n = 0; n < NumNodes; ++n)
		{
			int NumConnEdges = NodesCount[n];
			for (int k = 0; k < NumConnEdges; ++k)
			{
				int connNode = NtoNodes[n][k];
				if (connNode > n)
				{
					JunctionList[i]->EdgtoVert[EdgCount*2] = n;
					JunctionList[i]->EdgtoVert[EdgCount*2+1] = connNode;
					EdgCount += 1;
	
				}
		
			}
		
	
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
			JunctionList[i]->EdgtoEls[e*2] = -1;
			JunctionList[i]->EdgtoEls[e*2+1] = -1;
		}
		
		for (int e = 0; e < NumEdges; ++e)
		{
			int v0 = JunctionList[i]->EdgtoVert[e*2];
			int v1 = JunctionList[i]->EdgtoVert[e*2+1];
			int elfound = 0;
			for (int el = 0; el < NumEl; ++el)
			{
				for (int k = 0; k < 3; ++k)
				{
					int n0 = JunctionList[i]->EltoVert[el*3+k];
					int ind1 = (k+1)%3;
					int n1 = JunctionList[i]->EltoVert[el*3+ind1];
					if (((n0 == v0) && (n1 == v1))|| ((n0 == v1) &&(n1==v0)))
					{
						if (JunctionList[i]->EdgtoEls[e*2] < 0)
							JunctionList[i]->EdgtoEls[e*2] = el;
						else if (JunctionList[i]->EdgtoEls[e*2+1] < 0)
							JunctionList[i]->EdgtoEls[e*2+1] = el;
						elfound+=1;
						break;
					}
				}
				if (elfound==2)
					break;
				if ((el == (NumEl-1)) && elfound == 1)
					JunctionList[i]->EdgtoEls[e*2+1] = JunctionList[i]->EdgtoEls[e*2];
	
			 }
			if (elfound == 0)
			{
				printf("edge %d not an element edge \n",e);
				exit(EXIT_FAILURE);
			}
	
		}

		// Read and mark Inflow and Outflow edges and the channels associated with them	
		JunctionList[i]->BdryPrescribed = malloc(NumEdges*sizeof(int));
		JunctionList[i]->ChannelNumber = malloc(NumEdges*sizeof(int));

	
		double minEdgLength = 99999999999;
		for (int e = 0; e < NumEdges; ++e)
		{
			int n0 = JunctionList[i]->EdgtoVert[e*2+0];
			int n1 = JunctionList[i]->EdgtoVert[e*2+1];

			double edgLength = sqrt((JunctionList[i]->x[n1]-JunctionList[i]->x[n0])*(JunctionList[i]->x[n1]-JunctionList[i]->x[n0])+(JunctionList[i]->y[n1]-JunctionList[i]->y[n0])*(JunctionList[i]->y[n1]-JunctionList[i]->y[n0]));
			minEdgLength = fmin(minEdgLength,edgLength);

			int inflowOutflow = 0;
			for (int k = 0 ; k < NumInflowOutflow; ++k)
			{
				int b0 = InOutNodes[k][0];
				int b1 = InOutNodes[k][1];
				if (((b0==n0) && (b1 ==n1)) || ((b0==n1) && (b1==n0)))
				{
					JunctionList[i]->BdryPrescribed[e] = InOutNodes[k][3];
					JunctionList[i]->ChannelNumber[e] = InOutNodes[k][2];
					inflowOutflow = 1;
					break;
				}

			}

			if (!inflowOutflow)
			{
				JunctionList[i]->BdryPrescribed[e] = 0;	
				JunctionList[i]->ChannelNumber[e] = -1;
			}

		}
		JunctionList[i]->minEdgLength = minEdgLength;

	} // end loop junctions

	fclose(JuncNodes);

	// count the number of junction edges connected to each channel. 
	// Assuming that one channel cannot be connected to more than a hundred junction edges:
	int InflowConnList[NumChannels];
	int OutflowConnList[NumChannels];

	int InflowEdges[NumChannels][100];
	int OutflowEdges[NumChannels][100];

	int InflowJunctionNumber[NumChannels][100];
	int OutflowJunctionNumber[NumChannels][100];

	int InflowEdgNum[NumChannels];
	int OutflowEdgNum[NumChannels];

	for (int i = 0; i < NumChannels; ++i)
	{
		InflowConnList[i] = 0;
		OutflowConnList[i] = 0;
		InflowEdgNum[i] = 0;
		OutflowEdgNum[i] = 0;
		for (int j = 0; j < 100; ++j)
		{
			InflowEdges[i][j] = -1;
			OutflowEdges[i][j] = -1;
		}
	}

	for (int i = 0; i < NumJunctions; ++i)
	{
		int NumEdges = JunctionList[i]->TotalNumEdges;
		for (int j = 0; j < NumEdges; ++j)
		{
			int conn_to_channel = JunctionList[i]->BdryPrescribed[j];
			int channelNum = JunctionList[i]->ChannelNumber[j];

			
			// inflow for channel, outflow for junction
			if (conn_to_channel==2)
			{
				InflowConnList[channelNum] = InflowConnList[channelNum] + 1;
				InflowEdges[channelNum][InflowEdgNum[channelNum]] = j;
				InflowJunctionNumber[channelNum][InflowEdgNum[channelNum]] = i;
				InflowEdgNum[channelNum] = InflowEdgNum[channelNum] + 1;
			}
		
			else if (conn_to_channel == 1)
			{
				OutflowConnList[channelNum] = OutflowConnList[channelNum] + 1;
				OutflowEdges[channelNum][OutflowEdgNum[channelNum]] = j;
				OutflowJunctionNumber[channelNum][OutflowEdgNum[channelNum]] = i;
				OutflowEdgNum[channelNum] = OutflowEdgNum[channelNum] + 1;
			}
				
		}
	}

	for (int i =0; i < NumChannels; ++i)
	{
		ChannelList[i]->NumInflowJunctionEdges = InflowConnList[i];
		ChannelList[i]->NumOutflowJunctionEdges = OutflowConnList[i];

		ChannelList[i]->InflowJunctionEdges = malloc(2*InflowConnList[i]*sizeof(int));
		ChannelList[i]->OutflowJunctionEdges = malloc(2*OutflowConnList[i]*sizeof(int));

		int j = 0;
		while (InflowEdges[i][j] != -1)
		{
			ChannelList[i]->InflowJunctionEdges[2*j] = InflowEdges[i][j];
			ChannelList[i]->InflowJunctionEdges[2*j+1] = InflowJunctionNumber[i][j];
			++j;
		}

		int k = 0;
		while(OutflowEdges[i][k] != -1)
		{
			ChannelList[i]->OutflowJunctionEdges[2*k] = OutflowEdges[i][k];
			ChannelList[i]->OutflowJunctionEdges[2*k+1] = OutflowJunctionNumber[i][k];
			++k;
		}

	}

/*	FILE *juncCenter = fopen("JuncCenterNodes.in","r");
	int NumElx;
	char line1[500];
	fgets(line1, sizeof(line), juncCenter);
	sscanf(line1, "%d", &NumElx);
	
	FILE *writeFile = fopen("JuncCenterNodes.out","w");

	for (int i = 0; i < NumElx; ++i)
	{
		int node1, node2;
		char newline[500];
		fgets(newline, sizeof(newline), juncCenter);
		sscanf(newline, "%d\t%d", &node1, &node2);
		node1 = node1-1;
		node2 = node2-1;
		
		int el_found, edg;
		for (int el =0; el < JunctionList[0]->NumEl; ++el)
		{
			int v1 = JunctionList[0]->EltoVert[el*3];
			int v2 = JunctionList[0]->EltoVert[el*3+1];
			int v3 = JunctionList[0]->EltoVert[el*3+2];

			if ((node1 == v1 && node2 == v2) || (node1 == v2 && node2 == v1))
			{
				el_found = el;
				edg = 0;
				break;
			}
			else if ((node1 == v2 && node2 == v3) || (node1 == v3 && node2 == v2))
			{
				el_found = el;
				edg = 1;
				break;
			}
			else if ((node1 == v3 && node2 == v1) || (node1 == v1 && node2 == v3))
			{
				el_found = el;
				edg = 2;
				break;
			}
		}

		fprintf(writeFile, "%d\t%d\n", el_found, edg);
	}	

	fclose(writeFile);
	fclose(juncCenter);
*/			

}	


