#include <stdio.h>
#include <time.h>
#include <linux/limits.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "Globals.h"
#include "constitutive_equations.h"
#include "Watershed.h"
#include "manufactured_solution.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

extern void GetLGLWeights(int N, double* w);

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

static int firstOutput = 1;
static char folderName[500];
static char ChannelHeightFileName[PATH_MAX + 1];
static char ChannelQFileName[PATH_MAX + 1];
static char KinFileName[PATH_MAX + 1];
static char JuncHeightFileName[PATH_MAX + 1];
static char JuncQFileName[PATH_MAX + 1];
static char FpHeightFileName[PATH_MAX + 1];
static char FpQFileName[PATH_MAX + 1];

void outputFloodplainDataToFile(double currtime, double finalTime, double recordTimeIntervals)
{
	char openFileFormat[1];
	int writeHeader = 0;
	if (firstOutput)
	{
		time_t rawtime;
		time(&rawtime);
		struct tm *timeinfo = localtime(&rawtime);
		sprintf(folderName, "DataOutput_%d_%02d_%02d_%02d_%02d_%02d",timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
		if (mkdir(folderName, 0777) == -1)
			perror("The following error occured while making the data output directory");
		
		sprintf(FpHeightFileName, "%s/Floodplains.63", folderName);
		sprintf(FpQFileName, "%s/Floodplains.64", folderName);

		openFileFormat[0] = 'w';
		writeHeader = 1;
		firstOutput = 0;
	}
	else
	{
		openFileFormat[0] = 'a';
	}

	// Output Floodplain data
	FILE *file1 = fopen(FpHeightFileName, openFileFormat);
	FILE *file2 = fopen(FpQFileName, openFileFormat);
	if (!file1)
		printf("Error while opening file1\n");
	if (!file2)
		printf("Error while opening file2\n");
	if (writeHeader)
	{
		int numRecords;
		if (fmod(finalTime, recordTimeIntervals) == 0)
			numRecords = finalTime/recordTimeIntervals + 1;
		else
			numRecords = finalTime/recordTimeIntervals + 2;

		fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file1, "Number of Floodplains = %d\n", NumFloodplains);
		fprintf(file2, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file2, "Number of Floodplains = %d\n", NumFloodplains);
		for (int j = 0; j < NumFloodplains; j++)
		{
			fprintf(file1, "Number of data points for junction %d = %d\n", j, FloodplainList[j]->NumEl);
			fprintf(file2, "Number of data points for junction %d = %d\n", j, FloodplainList[j]->NumEl);
		}
	}
	for (int i=0; i < NumFloodplains; i++)
	{
		fprintf(file1, "# Average H and Zeta at time %.3f for Junction %d\n", currtime, i); 
		fprintf(file2, "# Average Qx and Qy at time %.3f for Junction %d\n", currtime, i);

		for (int j=0; j < FloodplainList[i]->NumEl; ++j)
		{
			double avgZeta = 0, avgHeight = 0, avgQx = 0, avgQy =0, avgz =0, zarr[3];
			for(int k=0; k<3; ++k)	
			{
				double zeta = FloodplainList[i]->zeta[j][k];
				double Qx = FloodplainList[i]->Qx[j][k];
				double Qy = FloodplainList[i]->Qy[j][k];
				zarr[k]= FloodplainList[i]->NodalZ[j][k];
				double height = zeta + zarr[k];
				
				avgZeta += zeta;
				avgHeight += height;
				avgQx += Qx;
				avgQy += Qy;
				avgz += zarr[k];
			}
			fprintf(file1, "%d \t %.13f \t %.13f \n", j, avgHeight/3, avgZeta/3);
			fprintf(file2, "%d \t %.13f \t %.13f\n", j, avgQx/3, avgQy/3 );
		}		
	}
	fclose(file1);
	fclose(file2);

}

void output2DNodalError(double time)
{
	char *fileName1 = "NodalErrorsH.out";
	char *fileName2 = "NodalErrorsQ.out";
	FILE* ef1 = fopen(fileName1, "w");
	FILE* ef2 = fopen(fileName2, "w");

	int i = 0;
	for (int j=0; j < FloodplainList[i]->NumVerts; ++j)
	{
		double xval = FloodplainList[i]->Vx[j];
		double yval = FloodplainList[i]->Vy[j];
		int num_conn_els = FloodplainList[i]->ElCount[i];
		double nodalZeta = 0;
		double nodalQx = 0;
		double nodalQy = 0;
		for (int k = 0; k < num_conn_els; k++)
		{
			int curr_el = FloodplainList[i]->VtoEl[j][k];
			int loc_v;
			// find the local vertex number for curr_el
			for (int v = 0; v < 3; v++)
			{
				if (j == FloodplainList[i]->EltoVert[curr_el*3+v])
				{
					loc_v = v;
					break;
				}
			}

			int loc_node = FloodplainList[i]->VtoNode[curr_el][loc_v];
			nodalZeta += FloodplainList[i]->zeta[curr_el][loc_node];
			nodalQx += FloodplainList[i]->Qx[curr_el][loc_node];
			nodalQy += FloodplainList[i]->Qy[curr_el][loc_node];

		}
		nodalZeta /= num_conn_els;
		nodalQx /= num_conn_els;
		nodalQy /= num_conn_els;

		double manzeta = getmanH(xval, yval, time);
		double manQx = getQx(xval, yval, time);
		double manQy = getQy(xval, yval, time);

		double er1 = -nodalZeta + manzeta;
		double er2 = -nodalQx + manQx;
		double er3 = -nodalQy + manQy;

		fprintf(ef1, "%d \t %.13f\n", j, er1);
		fprintf(ef2, "%d \t %.13f \t %.13f\n", j, er2, er3);
	}

	fclose(ef1);
	fclose(ef2);

}

void outputFloodplainNodalDataToFile(double currtime, double finalTime, double recordTimeIntervals)
{
	char openFileFormat[1];
	int writeHeader = 0;
	if (firstOutput)
	{
		time_t rawtime;
		time(&rawtime);
		struct tm *timeinfo = localtime(&rawtime);
		sprintf(folderName, "DataOutput_%d_%02d_%02d_%02d_%02d_%02d",timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
		if (mkdir(folderName, 0777) == -1)
			perror("The following error occured while making the data output directory");
		
		sprintf(FpHeightFileName, "%s/Floodplains.63", folderName);
		sprintf(FpQFileName, "%s/Floodplains.64", folderName);

		openFileFormat[0] = 'w';
		writeHeader = 1;
		firstOutput = 0;
	}
	else
	{
		openFileFormat[0] = 'a';
	}

	// Output Floodplain data
	FILE *file1 = fopen(FpHeightFileName, openFileFormat);
	FILE *file2 = fopen(FpQFileName, openFileFormat);
	if (!file1)
		printf("Error while opening file1\n");
	if (!file2)
		printf("Error while opening file2\n");
	if (writeHeader)
	{
		int numRecords;
		if (fmod(finalTime, recordTimeIntervals) == 0)
			numRecords = finalTime/recordTimeIntervals + 1;
		else
			numRecords = finalTime/recordTimeIntervals + 2;

		fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file1, "Number of Floodplains = %d\n", NumFloodplains);
		fprintf(file2, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file2, "Number of Floodplains = %d\n", NumFloodplains);
		for (int j = 0; j < NumFloodplains; j++)
		{
			fprintf(file1, "Number of data points for junction %d = %d\n", j, FloodplainList[j]->NumVerts);
			fprintf(file2, "Number of data points for junction %d = %d\n", j, FloodplainList[j]->NumVerts);
		}
	}
	for (int i=0; i < NumFloodplains; i++)
	{
		fprintf(file1, "# Nodal H at time %.3f for Junction %d\n", currtime, i); 
		fprintf(file2, "# Nodal Qx and Qy at time %.3f for Junction %d\n", currtime, i);

		for (int j=0; j < FloodplainList[i]->NumVerts; ++j)
		{
			int num_conn_els = FloodplainList[i]->ElCount[i];
			double nodalZeta = 0;
			double nodalQx = 0;
			double nodalQy = 0;
			for (int k = 0; k < num_conn_els; k++)
			{
				int curr_el = FloodplainList[i]->VtoEl[j][k];
				int loc_v;
				// find the local vertex number for curr_el
				for (int v = 0; v < 3; v++)
				{
					if (j == FloodplainList[i]->EltoVert[curr_el*3+v])
					{
						loc_v = v;
						break;
					}
				}

				int loc_node = FloodplainList[i]->VtoNode[curr_el][loc_v];
				nodalZeta += FloodplainList[i]->zeta[curr_el][loc_node];
				nodalQx += FloodplainList[i]->Qx[curr_el][loc_node];
				nodalQy += FloodplainList[i]->Qy[curr_el][loc_node];

			}
			nodalZeta /= num_conn_els;
			nodalQx /= num_conn_els;
			nodalQy /= num_conn_els;
			
			fprintf(file1, "%d \t %.13f\n", j, nodalZeta);
			fprintf(file2, "%d \t %.13f \t %.13f\n", j, nodalQx, nodalQy );
		}		
	}
	fclose(file1);
	fclose(file2);

}
void outputDataToFile(double currtime, double finalTime, double recordTimeIntervals)
{
	char openFileFormat[1];
	int writeHeader = 0;
	if (firstOutput)
	{
		time_t rawtime;
		time(&rawtime);
		struct tm *timeinfo = localtime(&rawtime);
		sprintf(folderName, "DataOutput_%d_%02d_%02d_%02d_%02d_%02d",timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
		if (mkdir(folderName, 0777) == -1)
			perror("The following error occured while making the data output directory");
		
		sprintf(ChannelHeightFileName, "%s/Channels.63", folderName);
		sprintf(ChannelQFileName, "%s/Channels.64", folderName);
		sprintf(KinFileName, "%s/KinField.63", folderName);
		sprintf(JuncHeightFileName, "%s/Junctions.63", folderName);
		sprintf(JuncQFileName, "%s/Junctions.64", folderName);
		//sprintf(FpHeightFileName, "%s/Floodplains.63", folderName);
		//sprintf(FpQFileName, "%s/Floodplains.64", folderName);

		// Move junction grid data file to this newly created folder
		/*char NewJunctionGridFileName[PATH_MAX + 1];
		sprintf(NewJunctionGridFileName, "%s/JunctionMesh.14", folderName);
		int ret = rename("./Output/JunctionMesh.14", NewJunctionGridFileName);
		if (ret != 0)
			printf("Error: unable to move the junction grid file\n");
*/
		// Move flowpaths file to this newly created folder
		int ret;
		char NewFlowPathsFileName[PATH_MAX + 1];
		sprintf(NewFlowPathsFileName, "%s/FlowPaths.out", folderName);
		ret = rename("./Output/FlowPaths.out", NewFlowPathsFileName);
		if (ret != 0)
			printf("Error: unable to move the flow paths file\n");

		// Move ChannelNodes.14 to this newly created folder
/*		char NewChannelsGridFileName[PATH_MAX + 1];
		sprintf(NewChannelsGridFileName, "%s/Channels.14", folderName);
		ret = rename("./Output/Channels.14", NewChannelsGridFileName);
		if (ret != 0)
			printf("Error: unable to move channels grid file\n");
*/
		openFileFormat[0] = 'w';
		writeHeader = 1;
		firstOutput = 0;
	}
	else
	{
		openFileFormat[0] = 'a';
	}

		//char cwd[PATH_MAX+1];
	//if (getcwd(cwd, sizeof(cwd))==NULL)
	//	perror("Error while getting current working directory");
	
	
	//Output Channels data 
	FILE* file1;
	FILE* file2;
	file1 = fopen(ChannelHeightFileName, openFileFormat);
	file2 = fopen(ChannelQFileName, openFileFormat);
	//if (!file1 || !file2)
	//{
	//	printf("Could not open file for writing. Exiting now\n");
	//	exit(1);
	//}
	if (writeHeader)
	{
		int numRecords;
		if (fmod(finalTime, recordTimeIntervals) == 0)
			numRecords = finalTime/recordTimeIntervals+1;
		else
			numRecords = finalTime/recordTimeIntervals + 2;
		fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file1, "Number of Channels = %d\n", NumChannels);
		fprintf(file2, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file2, "NumChannels = %d\n", NumChannels);
		for (int c = 0; c < NumChannels; c++)
		{
			fprintf(file1, "Number of data points for channel %d = %d\n", c, ChannelList[c]->NumNodes);
			fprintf(file2, "Number of data points for channel %d = %d\n", c, ChannelList[c]->NumNodes);
		}
	}
	for (int i=0; i<NumChannels; i++)
	{
		fprintf(file1, "# H and Zeta at time %.3f for Channel %d\n", currtime, i); 
		fprintf(file2, "#Q at time %.3f for Channel %d\n", currtime, i);

		int NumNodes = ChannelList[i]->NumNodes;	

		for (int j = 0; j<NumNodes; j++)
		{	
			double bval = ChannelList[i]->NodalB[j];
			double zval = ChannelList[i]->NodalZ[j];
			double x_val = ChannelList[i]->NodalX[j];
			double y_val = ChannelList[i]->NodalY[j];
			double A = ChannelList[i]->A[j+1]; 
			double Q = ChannelList[i]->Q[j+1];
			double m1val = ChannelList[i]->Nodalm1[j];
			double m2val = ChannelList[i]->Nodalm2[j];
			double H = getH(A, bval, m1val, m2val); 
			double zeta = H - zval;
			fprintf(file1, "%.13f \t %.13f\n", H, zeta);
			fprintf(file2, "%.13f\n", Q);

		}

	}

	fclose(file1);
	fclose(file2);

	// Output Junction data
	file1 = fopen(JuncHeightFileName, openFileFormat);
	file2 = fopen(JuncQFileName, openFileFormat);
	if (!file1)
		printf("Error while opening file1\n");
	if (!file2)
		printf("Error while opening file2\n");
	if (writeHeader)
	{
		int numRecords;
		if (fmod(finalTime, recordTimeIntervals) == 0)
			numRecords = finalTime/recordTimeIntervals + 1;
		else
			numRecords = finalTime/recordTimeIntervals + 2;

		fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file1, "Number of Junctions = %d\n", NumJunctions);
		fprintf(file2, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		fprintf(file2, "NumJunctions = %d\n", NumJunctions);
		for (int j = 0; j < NumJunctions; j++)
		{
			fprintf(file1, "Number of data points for junction %d = %d\n", j, JunctionList[j]->NumEl);
			fprintf(file2, "Number of data points for junction %d = %d\n", j, JunctionList[j]->NumEl);
		}
	}
	for (int i=0; i<NumJunctions; i++)
	{
		
		fprintf(file1, "# Average H and Zeta at time %.3f for Junction %d\n", currtime, i); 
		fprintf(file2, "# Average Qx and Qy at time %.3f for Junction %d\n", currtime, i);
		
		for (int j=0; j<JunctionList[i]->NumEl; ++j)
		{
			double avgZeta = 0, avgHeight = 0, avgQx = 0, avgQy =0, avgz =0, zarr[3];
			for(int k=0; k<3; ++k)	
			{
				double zeta = JunctionList[i]->zeta[j][k];
				double Qx = JunctionList[i]->Qx[j][k];
				double Qy = JunctionList[i]->Qy[j][k];
				zarr[k]= JunctionList[i]->NodalZ[j][k];
				double height = zeta + zarr[k];
				
				avgZeta += zeta;
				avgHeight += height;
				avgQx += Qx;
				avgQy += Qy;
				avgz += zarr[k];

			}
			fprintf(file1, "%d \t %.13f \t %.13f \n", j, avgHeight/3, avgZeta/3);
			fprintf(file2, "%d \t %.13f \t %.13f\n", j, avgQx/3, avgQy/3 );
		}		
	}
	fclose(file1);
	fclose(file2);

	// Output Kinematic Field Data
	// Here we assume that we are working with first order polynomials only
	if (FloodplainList[0]->floodedStatus == 0)
	{
		file1 = fopen(KinFileName, openFileFormat);
		if (writeHeader)
		{
			int numRecords;
			if (fmod(finalTime, recordTimeIntervals) == 0)
				numRecords = finalTime/recordTimeIntervals + 1;
			else
				numRecords = finalTime/recordTimeIntervals + 2;

			fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
		}
		fprintf(file1,"# Average height at time %3.3f in Kinematic Field\n", currtime);

		int NumEl = FloodplainList[0]->NumEl;
		for (int i = 0; i < NumEl; i++)
		{
			if (KinematicElList[i]->isActive == 1)
			{
				double weq = KinematicElList[i]->weq;
				double height = KinematicElList[i]->A[0]/weq;
				if (KinematicElList[i]->numUpstreamEls > 0)
				{
					int numUpstreamEls = KinematicElList[i]->numUpstreamEls;
					height = KinematicElList[i]->A[0]/weq;
					for (int j = 0; j < numUpstreamEls; j++)
					{
						int el = KinematicElList[i]->upstreamEls[j];
						height += KinematicElList[el]->A[1]/weq;
					}
					height = height/(numUpstreamEls+1);
				}
				//printf("%d %.13f\n", i, height);
				fprintf(file1, "%d %.13f\n", i, height);
			}
			else
				fprintf(file1, "%d %1.13f\n", i, 1e-7);
		}
		fclose(file1);
	}

	//static int previouslyFlooded = 0;
	//static char fpOpenFileFormat[1];
	//int writefpHeader = 0;

	//// Output Floodplain data if the channels have flooded
	//if (FloodplainList[0]->floodedStatus == 1)
	//{
	//	if (previouslyFlooded)
	//	{
	//		fpOpenFileFormat[0] = 'a';
	//	}
	//	else
	//	{
	//		fpOpenFileFormat[0] = 'w';
	//		writefpHeader = 1;	
	//	}

	//	file1 = fopen(FpHeightFileName, fpOpenFileFormat);
	//	file2 = fopen(FpQFileName, fpOpenFileFormat);

	//	
	//	if (!file1)
	//		printf("Error while opening file1\n");
	//	if (!file2)
	//		printf("Error while opening file2\n");
	//	if (writefpHeader)
	//	{
	//		int numRecords;
	//		if (fmod(finalTime, recordTimeIntervals) == 0)
	//			numRecords = finalTime/recordTimeIntervals + 1;
	//		else
	//			numRecords = finalTime/recordTimeIntervals + 2;

	//		fprintf(file1, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
	//		fprintf(file2, "End of Simulation Time = %lf.\nTotal Number of Records = %d\n", finalTime, numRecords);
	//		fprintf(file1, "Number of data points for floodplain = %d\n", FloodplainList[0]->NumEl);
	//		fprintf(file2, "Number of data points for floodplain = %d\n", FloodplainList[0]->NumEl);
	//	}
	//		
	//	fprintf(file1, "# Average H and Zeta at time %.3f\n", currtime); 
	//	fprintf(file2, "# Average Qx and Qy at time %.3f\n", currtime);
	//
	//	for (int j=0; j<FloodplainList[0]->NumEl; ++j)
	//	{
	//		double avgZeta = 0, avgHeight = 0, avgQx = 0, avgQy =0, avgz =0, zarr[3];
	//		for(int k=0; k<3; ++k)	
	//		{
	//			double zeta = FloodplainList[0]->zeta[j][k];
	//			double Qx = FloodplainList[0]->Qx[j][k];
	//			double Qy = FloodplainList[0]->Qy[j][k];
	//			zarr[k]= FloodplainList[0]->NodalZ[j][k];
	//			double height = zeta + zarr[k];
	//			
	//			avgZeta += zeta;
	//			avgHeight += height;
	//			avgQx += Qx;
	//			avgQy += Qy;
	//			avgz += zarr[k];

	//		}
	//		fprintf(file1, "%d \t %.13f \t %.13f \n", j, avgHeight/3, avgZeta/3);
	//		fprintf(file2, "%d \t %.13f \t %.13f\n", j, avgQx/3, avgQy/3 );
	//	}		
	//	fclose(file1);
	//	fclose(file2);

	//	previouslyFlooded = 1;

	//}

}

void calculate_l2_error_2d(double time)
{
	int NumEl = FloodplainList[0]->NumEl;
	int Np = FloodplainList[0]->Np;
	double errsqH = 0, errsqQx = 0, errsqQy = 0;

	gsl_vector *diffH = gsl_vector_alloc(Np);
	gsl_vector *diffQx = gsl_vector_alloc(Np);
	gsl_vector *diffQy = gsl_vector_alloc(Np);

	//double max = 0;
	for (int i = 0; i < NumEl; i++)
	{
		for (int j = 0; j < Np; j++)
		{
			double xval = FloodplainList[0]->NodalX[i][j];
			double yval = FloodplainList[0]->NodalY[i][j];
			double exactH = getmanH(xval, yval, time);
			double exactQx = getQx(xval, yval, time);
			double exactQy = getQy(xval, yval, time);
			double compH = FloodplainList[0]->zeta[i][j];
			double compQx =  FloodplainList[0]->Qx[i][j];
			double compQy = FloodplainList[0]->Qy[i][j];
			gsl_vector_set(diffH, j, fabs(exactH-compH));
			gsl_vector_set(diffQx, j, fabs(exactQx-compQx));
			gsl_vector_set(diffQy, j, fabs(exactQy-compQy));
			//gsl_vector_set(diffH, j, RHSZeta[i*Np+j]);
			//gsl_vector_set(diffQx, j,  RHSQx[i*Np+j]);
			//gsl_vector_set(diffQy,j, RHSQy[i*Np+j]);
			//max = fmax(fabs(exactQx - compQx), max);
					
		}

		gsl_vector *tmpH = gsl_vector_alloc(Np);
		gsl_vector *tmpQx = gsl_vector_alloc(Np);
		gsl_vector *tmpQy = gsl_vector_alloc(Np);
		double jac = FloodplainList[0]->jac[i];
		gsl_blas_dgemv(CblasNoTrans, jac, MassMatrix2D, diffH, 0.0, tmpH);
		gsl_blas_dgemv(CblasNoTrans, jac, MassMatrix2D, diffQx, 0.0, tmpQx);
		gsl_blas_dgemv(CblasNoTrans, jac, MassMatrix2D, diffQy, 0.0, tmpQy);

		double elerrHsq, elerrQxsq, elerrQysq;
		gsl_blas_ddot(diffH, tmpH, &elerrHsq);
		gsl_blas_ddot(diffQx, tmpQx, &elerrQxsq);
		gsl_blas_ddot(diffQy, tmpQy, &elerrQysq);

		errsqH += elerrHsq;
		errsqQx += elerrQxsq;
		errsqQy += elerrQysq;

		gsl_vector_free(tmpH);
		gsl_vector_free(tmpQx);
		gsl_vector_free(tmpQy);
	}

	gsl_vector_free(diffH);
	gsl_vector_free(diffQx);
	gsl_vector_free(diffQy);

	printf("time at error computation = %lf\n", time);
	printf("l2 error in H = %lf \n", sqrt(errsqH));
	printf("l2 error in Qx = %lf \n", sqrt(errsqQx));
	printf("l2 error in Qy = %lf \n", sqrt(errsqQy));
	//printf("maximum difference in Qx = %lf\n", max);
	
}

void append_to_file(char *fileName, double time, double data)
{
	FILE* file = fopen(fileName, "a");
	fprintf(file, "%lf %lf\n", time,data);
	fclose(file);
}


double calculateTotalWater(struct channel *Chan)
{
	int NumEl = Chan->NumEl;
	int Np = Chan->Np;
	int P = Chan->P;

	double totalWater = 0;
	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);
	for (int k = 0; k < NumEl; k++)
	{
		double water = 0;
		double dh = Chan->dh[k];
		for (int i = 0; i < Np; i++)
		{
			double m1val = Chan->Nodalm1[k*Np+i];
			double m2val = Chan->Nodalm2[k*Np+i];
			double bval = Chan->NodalB[k*Np+i];
			double A = Chan->A[k*Np+i+1];
			double h = getH(A, bval, m1val, m2val);
			water += LGLWeight[i]*(h-0.0000001);
		}
		totalWater += 0.5*dh*water;
	}

	return totalWater;

}


double ftTom(double ft)
{
	double meter = ft*0.3048;
	return meter;
}

double mToft (double m)
{
	double ft = m*3.28084;
	return ft;

}

