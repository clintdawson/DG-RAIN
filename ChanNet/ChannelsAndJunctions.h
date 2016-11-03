
/*****************************************************************************//**
 * @file ChannelsAndJunctions.h
 * This file contains channel and junction structure definitions.               
 ******************************************************************************/                                      
 
/***************************************************************************//**
 * A structure used to store all the variables associated with
 *  the discretization and the physical properties of a channel. Other than the 
 *  solutions (A and Q), the maximum eigenvalue (max_lambda) and the wet/dry status
 *  (WD), every other field is either obtained or can be determined from the grid 
 *  file.                   
 *******************************************************************************/  
 #ifndef CHANNELS_JUNCTIONS

#define CHANNELS_JUNCTIONS

struct channel 
{
	int NumEl;          ///< Number of elements used for the discretization of the channel
						///< 
	int NumNodes;       ///< Total number of nodes (counts one node only once)
						///<
	int* GlobalNodeNum; ///< Global node numbers of the nodes of the channel
						 ///<
	
	/* fields to store the physical attributes of the channel*/
	double *x;          ///< x-coordinates of the nodes of the channel
						///<
	double *y;          ///< y-coordinates of the nodes of the channel
						///<
	double *b;          ///< width of the channel at each node
						///<
	double *z;          ///< bathymetry at each node
						///<
	
	double *nFriction;  ///< manning's coefficient for the channel prescribed at each node
						///<
	
//	int JunctionNode[2];    ///< store the global node number of the node connected  to a junction
	
	double mindh;       ///< size of the smallest element of the channel.
						///<
	double *beta;       ///< The momentum correction coefficient for the channel prescribed
	                    ///< at each node
	
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
	                             
	
	/* fields to store the solution */
	double *A;          ///< Stores the wet cross-sectional area solution. 
	                    ///< Contains the solution at each node. Size = 2*NumNodes
	double *Q;          ///< Stores the volumetric discharge solution.
	                    ///< Contains the solution at each node. Size = 2*NumNodes
	double max_lambda;  ///< Maximum eigenvalue of the jacobian across all elements.
	                    ///< This information is recalculated at each time step and is
	                    ///< used to set the next time step according to the CFL condition
	
	
	/* Wetting and drying treatment */
	int *WD;            ///< Stores the wet/dry status of each element. Size = NumNodes-1
						///<
	        
};

/**********************************************************************************************//**
* A structure used to store all the variables associated with the discretization and the physical 
* properties of a 2-D junction element. Other than the solutions (zeta, Qx and Qy), the maximum
* eigenvalue (max_lambda) and the wet/dry status (WD), every other field is either obtained or can
* be determined from the grid file.
 ***********************************************************************************************/
struct junction
{
	int TotalNumEdges;              ///< Total number of the edges of the junction element.
	                                ///< Each edge is counted only once.
	int NumEl;                      ///< Number of elements in the junction
									///<
	int NumNodes;                   ///< Number of nodes in the junction. Each node is counted
	                                ///< only once.
	
	double *x;                      ///< x-coordinates of the nodes of the junction 
									///<
	double *y;                      ///< y-coordinates of the nodes of the junction
									///<
	double *z;                      ///< bathymetry at the nodes
									///<
	
	double minEdgLength;            ///< length of the smallest edge of an element.
									///<
	
	int *EdgtoEls;                  ///< In its row, EdgtoEls stores the elements that are 
	                                ///< connected by an edge. Size = TotalNumEdges X 2; 
	                                ///< (i,0) and (i,1) elements are the two elements connected
	                                ///< by edge i. If edge i is a boundary edge then it stores
	                                ///< the same element in both columns. 
	
	int *EltoVert;                  ///< Stores the global vertex number of the vertices of
	                                ///< the elements. size = NumEl X 3. (i,j)th element is the
	                                ///< global vertex number of the jth vertex of element i.
	
	int *EdgtoVert;                 ///< Stores the two vertices connected by the edges.
	                                ///< size = totalNumEdges X 2; (i,0) and (i,1) are the global
	                                ///< vertex numbers of the vertices connected by edge i.
	
	int *BdryPrescribed;            ///< For each edge, stores the information about whether 
	                                ///< or not the edge is connected to a channel. 
	                                ///< size = TotalNumEdges x 1. BdryPrescribed(i) is 1 if edge
	                                ///< i is connected to an inflow channel, 2 if edge i is
	                                ///< connected to an outflow channel and 0 if its neither.
	                                
	int *ChannelNumber;             ///< Stores the channel number associated with each edge. 
	                                ///< size = TotalNumEdges x 1. 
	
	//int AssocNode;
	//int AssocNodeLocalNumber;     // local node number of the node associated with the junction
	
	double *zeta;               ///< Stores the water surface height solution. size = NumEl x 3;
								///<
	double *Qx;                 ///< Stores the momentum in the x-direction. size = NumEl X 3;
								///<
	double *Qy;                 ///< Stores the momentum in the y-direction. size = NumEl X 3;
								///<
	
	double *bzeta;              ///< Stores the value of zeta prescribed at an edge. 
	                            ///< size = TotalNumEdges x 1. 0 if not connected to a channel.
	double *bQn;                ///< Stores the value of normal flow prescribed at an edge.
	                            ///< size = TotalNumEdges x 1. 9999 if not connected to a channel. 
	
	double *nFriction;          ///< Stores the Manning's coefficient at each node. size = NumNodes
								///<
	int *WD;                    ///< Stores the wet/dry status of each element. size = NumEl
								///<
	
	double max_lambda;          ///< Maximum eigenvalue of the jacobian across all elements.
	                            ///< This information is recalculated at each time step and is
	                            ///< used to set the next time step according to the CFL condition
};
    
#endif

