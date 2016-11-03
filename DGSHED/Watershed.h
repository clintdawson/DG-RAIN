#ifndef WATERSHED

#define WATERSHED

struct kinematicEl
{
  int el1, el2;
  double x1, x2;
  double y1, y2;
  double z1, z2;
	double nf1, nf2;
	double dh;			// length of the kinematic element
	double weq;			// equivalent width

  double node1Area;
  double node2Area;

  int connChanNum;
  int connChanEl;
  int isActive;
	
  int numUpstreamEls;
  int numDownstreamEls;
  int upstreamEls[2];

	int P;				// even though there is an option to use different P for different kinematic element,
								// the current implementation assumes that the P for each kinematic element 
								// is the same.
								// This is important in determining where in the RHS array values are stored. 
								// But this can be easily changed later
	int Np;
	double* NodalX;
	double* NodalY;
	double* NodalZ;
	double* NodalnFriction;
	double* dz;

	double* A;
	
};


struct channel
{
	int NumEl;
	int NumEdges;
	int P;
	int Np;
	int NumNodes;

	double* x;
	double* y;
	double* z;
	double* b;
	double* depth;
	double* m1;
	double* m2;

	double* dh;
	double* nFriction;

	double* NodalX;
	double* NodalY;
	double* NodalZ;
	double* NodalB;
	double* NodalDepth;
	double* Nodalm1;
	double* Nodalm2;
	double* NodalnFriction;

	double* db;
	double* dm1;
	double* dm2;
	double* dz;

	double* A;
	double* Q;
	double* qL;
	double* qM;
	
	double dt;
	double max_lambda;
	double mindh;

	int* originalGlobalNodeNum;  // only to be used to store connectivity between kinematic elements and channel elements.
															 // some node numbers are not accurate because they are modified while creating junctions

	double* beta;

	int connectedFlowFieldEl;

	double channelLength;

	int NumInflowJunctionEdges;  ///< Number of junction edges connected to the inflow end
	                             ///< (represented by the first node) of the channel
	int *InflowJunctionEdges;    ///< NumInflowJunctionEdges X 2 array stored as a 1-D array
	                             ///< in row-major order.
	                             ///< For row i, element (i,0) stores the (local) junction
	                             ///< edge number connected to the inflow end 
	                             ///<(represented by the first node) of the channel.
	                             ///< Element (i,1) stores the junction number connected to
	                             ///< the channel. 
	                
	int NumOutflowJunctionEdges; ///< Number of junction edges connected to the outflow end
	                             ///< (represented by the last node) of the channel 
	int *OutflowJunctionEdges;   ///< NumOutflowJunctionEdges X 2 array stored as a 1-D array
	                             ///< in row-major order.
	                             ///< For row i, element (i,0) stores the (local) junction edge
	                             ///< number connected to the outflow end (represented by the 
	                             ///< last node) of the channel.
	                             ///< Element (i,1) stores the junction number connected to 
	                             ///< the channel.
	
	/* Wetting and drying treatment */
	int *WD;            ///< Stores the wet/dry status of each element. Size = NumNodes-1
						///<
	
	int *BCtype; 				//< type of boundary conditions on the two boundaries of the channel. (size = 2)
	double **BCVals;

	int *NfloodBanks;   ///< stores if there is a bank that can be flooded on one or both sides of the channel. size = NumEl
	int *FloodPlainElNum; ///< stores which floodplain element a particular element of the channel is connected to

};

/**********************************************************************************************//**
* A structure used to store all the variables associated with the discretization and the physical 
* properties of a 2-D junction element. Other than the solutions (zeta, Qx and Qy), the maximum
* eigenvalue (max_lambda) and the wet/dry status (WD), every other field is either obtained or can
* be determined from the grid file.
 ***********************************************************************************************/
struct TwoDRegion
{
	int type;												////< Either 1 for junction or 2 for Floodplains
	int P;													/// < Order of polynomial approximation
	int Np;													/// < Number of nodes per element
	int TotalNumEdges;              ///< Total number of the edges of the junction element.
	                                ///< Each edge is counted only once.
	int NumEl;                      ///< Number of elements in the junction
									///<
	int NumVerts;                   ///< Number of vertices in the junction. 
	
	double *Vx;                      ///< x-coordinates of the vertices of the junction 
									///<
	double *Vy;                      ///< y-coordinates of the Vertices of the junction
									///<
	double *Vz;                      ///< bathymetry at the vertices
									///<
	
	double *VnFriction;          ///< Stores the Manning's coefficient at each vertex

	double **NodalX;						///< Stores the x-coordinates at the nodes in each element. size = NumEl X Np
	double **NodalY;						///< Stores the y-coordinates at the nodes in each element.
	double **NodalZ;						///< Stores the z-values at the nodes in each element.
	double **Nodaldzx;						///< Stores the dz/dx-values at the nodes in each element.
	double **Nodaldzy;						///< Stores the dz/dy-values at the nodes in each element.
	double **NodalnFriction; 		///< Stores the Manning's n values at the nodes in each element.
	
	double *jac; 								///< Stores the element jacobian
	double *edgJac;							///< Stores the jacobian for each edge. size = NumEl x 3
	double *nx;									///< size = NumEl*3; nx[3*i+j] is the x-component of the outward unit normal vector to edge j of element i
	double *ny;									///< size = NumEl*3; defined similarly to nx
	double *rx;									///< size = NumEl; For each element, store dr/dx.
	double *sx;									///< size = NumEl; For each element, store ds/dx.
	double *ry;									///< size = NumEl; For each element, store dr/dy.
	double *sy;									///< size = NumEl; For each element, store ds/dy.

	double minEdgLength;            ///< length of the smallest edge of an element.
									///<
	
	int *EdgtoEls;                  ///< In its row, EdgtoEls stores the elements that are 
	                                ///< connected by an edge. Size = TotalNumEdges X 2; 
	                                ///< (i,0) and (i,1) elements are the two elements connected
	                                ///< by edge i. If edge i is a boundary edge then it stores
	                                ///< the same element in both columns. 
	int *GlobaltoLocalEdg;					///< size = NumEdges X 2. For each edge, stores the local edge number for the two element it separates. 
	int **GlobalEdgPosNegNodes;			///< size = NumEdges X (Nfp x 2) . For each local node number on one side of the edge, store the local node number
																	///< for that same node on the opposite side of that edge

	int **PosInFVec;
	int *EltoVert;                  ///< Stores the global vertex number of the vertices of
	                                ///< the elements. size = NumEl X 3. (i,j)th element is the
	                                ///< global vertex number of the jth vertex of element i.
	int *EltoEdg;										///< This is only used to create kinematic elements from floodplains

	int *EdgtoVert;                 ///< Stores the two vertices connected by the edges.
	                                ///< size = totalNumEdges X 2; (i,0) and (i,1) are the global
	                                ///< vertex numbers of the vertices connected by edge i.
	
	int *BdryPrescribed;            ///< For each edge, stores the information about whether 
	                                ///< or not the edge is connected to a channel. 
	                                ///< size = TotalNumEdges x 1. BdryPrescribed(i) is 1 if edge
	                                ///< i is connected to an inflow channel, 2 if edge i is
	                                ///< connected to an outflow channel and 0 if its neither.
	                                
	int *ElCount;										///< size = NumVerts; Stores the number of elements sharing a vertex
	int **VtoEl;										///< In row i, stores the element that share vertex i
	int **VtoNode;								///<  size = NumEl x 3. In location [i,j], stores the local node number of vertex j of element i
	int *ChannelNumber;             ///< Stores the channel number associated with each edge. 
	                                ///< size = TotalNumEdges x 1. 
	int *ChannelElNumber;						///< Stores the channel element connected to each edge
	
	//int AssocNode;
	//int AssocNodeLocalNumber;     // local node number of the node associated with the junction
	
	double **zeta;               ///< Stores the water surface height solution. size = NumEl x Np;
								///<
	double **Qx;                 ///< Stores the momentum in the x-direction. size = NumEl X Np;
								///<
	double **Qy;                 ///< Stores the momentum in the y-direction. size = NumEl X Np ;
								///<
	
	double *bzeta;              ///< Stores the value of zeta prescribed at an edge. 
	                            ///< size = TotalNumEdges x 1. 0 if not connected to a channel.
	double *bQn;                ///< Stores the value of normal flow prescribed at an edge.
	                            ///< size = TotalNumEdges x 1. 9999 if not connected to a channel. 
	
								///<
	int *WD;                    ///< Stores the wet/dry status of each element. size = NumEl
								///<
	
	double max_lambda;          ///< Maximum eigenvalue of the jacobian across all elements.
	                            ///< This information is recalculated at each time step and is
	                            ///< used to set the next time step according to the CFL condition
															
	int floodedStatus;					///< Only used in floodplains. 1 if the connected channels have flooded. 0 if the connected channels have not flooded.
	int nChEdgs;									///< Stores the number of edges connected to channels
};

#endif
