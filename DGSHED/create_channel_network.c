//#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Watershed.h"
#include "Globals.h"
#include "process_twoD_domains.h"

/*********************************************************************************************************//**
* @file create_channel_network.c
* 
* This file contains code to create junctions from the channel nodes files and the full mesh file.
* The write_junction_mesh function interfaces with python for creating and meshing junctino.
* Then, channel networks are created from the grid files provided for the channels and the 
* junctions. 
*
* *********************************************************************************************************/
int NumChannels;		// Total number of channels
struct channel **ChannelList;

extern void *xcalloc(int items, int size);

void create_channel_network(char *ChannelCoordinates, char *JunctionCoordinates)
{
	FILE *ChanCoordinates = fopen(ChannelCoordinates, "r");
	char line[100];
	fgets(line, sizeof(line),ChanCoordinates);	
	sscanf(line, "%d", &NumChannels);

	// allocate space for an array to store channels
	ChannelList = malloc(NumChannels*sizeof(struct channel*));

	for (int i = 0; i < NumChannels; ++i)
	{
		char newline[100];
		int NumNodes;
		fgets(newline, sizeof(newline), ChanCoordinates);
		sscanf(newline, "%d", &NumNodes);

		double x[NumNodes], y[NumNodes], z[NumNodes], b[NumNodes], n[NumNodes], m1[NumNodes], m2[NumNodes];
	
		for (int j = 0; j < NumNodes; ++j)
		{
			char grid[1500];
			fgets(grid, sizeof(grid), ChanCoordinates);
			sscanf(grid, "%lf %lf %lf %lf %lf %lf %lf", &x[j], &y[j], &z[j], &b[j], &m1[j], &m2[j], &n[j]);

		}

		ChannelList[i] = malloc(sizeof(struct channel));

		ChannelList[i]->NumEdges = NumNodes;
		ChannelList[i]->NumEl = NumNodes-1;

		ChannelList[i]->x = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->y = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->z = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->depth = malloc(NumNodes*sizeof(double));
		ChannelList[i]->b = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->m1 = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->m2 = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->nFriction = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->dh = xcalloc(NumNodes-1,sizeof(double));
		
		for (int j = 0; j < NumNodes; ++j)
		{
			ChannelList[i]->x[j] = x[j];
			ChannelList[i]->y[j] = y[j];
			ChannelList[i]->z[j] = z[j];
			ChannelList[i]->b[j] = b[j];
			ChannelList[i]->m1[j] = m1[j];
			ChannelList[i]->m2[j] = m2[j];
			ChannelList[i]->nFriction[j] = n[j];
			ChannelList[i]->depth[j] = 10000000;		// use a really large number at first to
																							// indicate that the channel can hold an
																							// infinite amount of water
		}

		double chanLength = 0.0;
		double mindh = 9999999999;
		for (int j = 0; j < NumNodes-1; ++j)
		{
			double dh = sqrt((x[j+1]-x[j])*(x[j+1]-x[j])+(y[j+1]-y[j])*(y[j+1]-y[j]));
			ChannelList[i]->dh[j] = dh;
			mindh = fmin(mindh, dh);
			chanLength += dh;
		}	
	
		ChannelList[i]->mindh = mindh;
		ChannelList[i]->channelLength = chanLength;

	} // end loop for Channels

	// store channel boundary information
	// throw away first line
	char taline[100];
	fgets(taline, sizeof(taline), ChanCoordinates);
	for (int i = 0; i < NumChannels; i++)
	{
		ChannelList[i]->BCtype = xcalloc(2,sizeof(int));
		ChannelList[i]->BCVals = malloc(2*sizeof(int*));
		ChannelList[i]->BCVals[0] = xcalloc(2,sizeof(int));
		ChannelList[i]->BCVals[1] = xcalloc(2,sizeof(int));
	
	// initialize BCtype to -1
		ChannelList[i]->BCtype[0] = -1;
		ChannelList[i]->BCtype[1] = -1;

		int bc0, bc1;
		char g1[100];
		fgets(g1, sizeof(g1), ChanCoordinates);
		sscanf(g1, "%d %d", &bc0, &bc1);
		ChannelList[i]->BCtype[0] = bc0;
		ChannelList[i]->BCtype[1] = bc1;
		
		// if only height or discharge provided at the upstream end
		if (bc0 == 1 || bc0 ==2)
		{
			char g2[500];
			fgets(g2, sizeof(g2), ChanCoordinates);
			double val;
			sscanf(g2, "%lf", &val);
			if (bc0 ==1)
				ChannelList[i]->BCVals[0][0] = val;
			else
				ChannelList[i]->BCVals[0][1] = val;
		}
		
		// if both height and discharge provided at the upstream end
		else if (bc0 == 3)
		{
			char g2[500];
			fgets(g2, sizeof(g2), ChanCoordinates);
			double val1, val2;
			sscanf(g2, "%lf %lf", &val1, &val2);
			ChannelList[i]->BCVals[0][0] = val1;
			ChannelList[i]->BCVals[0][1] = val2;
		}
		
		// bctype 0 means no flow, btype 100 means connected to bay
		// btype 400 means connected to junction. 
		else if ((bc0 != 0) && (bc0 != 100) && (bc0 != 400) && (bc0 != 4))
		{
			printf("Unknown upstream boundary type for channel %d. Quitting the program now.\n", i);
			exit(1);
		}

		// if only height or discharge provided at the downstream end
		if (bc1 == 1 || bc1 ==2)
		{
			char g2[500];
			fgets(g2, sizeof(g2), ChanCoordinates);
			double val;
			sscanf(g2, "%lf", &val);
			if (bc0 ==1)
				ChannelList[i]->BCVals[1][0] = val;
			else
				ChannelList[i]->BCVals[1][0] = val;
		}
		
		// if both height and discharge provided at the downstream end
		else if (bc1 == 3)
		{
			char g2[500];
			fgets(g2, sizeof(g2), ChanCoordinates);
			double val1, val2;
			sscanf(g2, "%lf %lf", &val1, &val2);
			ChannelList[i]->BCVals[1][0] = val1;
			ChannelList[i]->BCVals[1][0] = val2;
		}
		
		// bctype 0 means no flow, btype 100 means connected to bay
		// btype 400 means connected to junction. 
		else if ((bc1 != 0) && (bc1 != 100) && (bc1 != 400) && (bc1 != 4))
		{
			printf("Unknown downstream boundary type for channel %d. Quitting the program now.\n", i);
			exit(1);
		}

	}

	fclose(ChanCoordinates);

	// now check to see if there are still channels with unspecified boundary conditions
	// on either edge
	for (int k = 0; k < NumChannels; k++)
	{
		if (ChannelList[k]->BCtype[0] == -1)
		{
			printf("Unspecified boundary condition on the first node of channel %d\n", k);
			printf("Imposing no flow condition\n");
			ChannelList[k]->BCtype[0] = 0;
		}

		if (ChannelList[k]->BCtype[1] == -1)
		{
			printf("Unspecified boundary condition on the last node of channel %d\n", k);
			printf("Imposing no flow condition\n");
			ChannelList[k]->BCtype[1] = 0;
		}
	}

		// read and process the grid file
	read2DGridFile(JunctionCoordinates, "Junction");
	
	for (int i = 0; i < NumJunctions; i++)
	{
		int NumEdges = JunctionList[i]->TotalNumEdges;

		// Create EdgtoEls
		JunctionList[i]->EdgtoEls = malloc(2*NumEdges*sizeof(int));
		createEdgtoEls(JunctionList[i]->EltoVert, JunctionList[i]->EdgtoVert, JunctionList[i]->NumEl, JunctionList[i]->TotalNumEdges, JunctionList[i]->EdgtoEls);
		
		// Create GlobaltoLocalEdg
		JunctionList[i]->GlobaltoLocalEdg = xcalloc(2*NumEdges, sizeof(int));
		createGlobaltoLocalEdg(JunctionList[i]->EltoVert, JunctionList[i]->EdgtoVert, JunctionList[i]->EdgtoEls, JunctionList[i]->TotalNumEdges, JunctionList[i]->GlobaltoLocalEdg);
		
		// create VToEl
		JunctionList[i]->VtoEl = malloc(JunctionList[i]->NumVerts*sizeof(int*));
		createVtoEl(JunctionList[i]->EltoVert, JunctionList[i]->ElCount, JunctionList[i]->NumEl, JunctionList[i]->NumVerts, JunctionList[i]->VtoEl);		
		
	} // end loop junctions


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
		
		// 400 means connected to junction
	/*	if (InflowConnList[i] > 0)
			ChannelList[i]->BCtype[0] = 400;
		if (OutflowConnList[i] > 0)
			ChannelList[i]->BCtype[1] = 400;
*/

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

}

