#pragma once
#include<iostream>
#include<fstream>
#include<algorithm>
#include<numeric>
#include<vector>
#include<string>
#include<chrono>
#include<ctime>
#include<cmath>
#include<limits>
#include<cstdio>
#include<queue>
#include <stack> 
#include <list>
#include <set>
#include<unordered_set> 
#include <iomanip>
#include<unordered_map>
#include<sstream>
using namespace std;

struct Node {
	double elevation;
	int nodeIndex; // original index from raster
	vector<int> parentsID;  ///Contour Tree parents node indices
	vector<int> childrenID;  //Contour Tree child

	// TODO: added by Saugat for connecting nearby trees
	vector<int> connectingTreeId;
	double curGain;
							 //contour tree construction
	vector<int> joinParentsID;
	int joinChildID;
	int splitParentID;
	vector<int> splitChildrenID;
	int upDegree; // for tree merge convinence
	int downDegree;
	int bank;
	double cost;
	float p;
	int originalId;
	int isObserved=0;
	int roughness;
	//collapse tree
	vector<int> collapsedPixelIDs;
	int isNa = 0;
	int isPits = 0;
	int isTree = 0;
	double fel;
	
	//message propagation
	//leaves to root
	vector<double> fi_ChildList;
	double fi[cNum];
	double fo[cNum];
	int foNode;
	int foNode_ischild;

	//root to leaves
	vector<double> go_ChildList;
	double go_parent[cNum];
	double gi[cNum];
	bool go_fromParent;
	bool go_fromChild;
	bool go_fromParentofChild;
	int go_lastVisitedNode;
	bool visited = false;
	//inference
	//added by Arpan
	vector<int> correspondingNeighbour;
	vector<int> correspondingNeighbourClassOne;
	vector<int> correspondingNeighbourClassZero;
	Node() {
	}
	Node(double value, int nodeindex) {
		elevation = value;
		nodeIndex = nodeindex;
		splitParentID = -1;
		joinChildID = -1;

	}
};



struct sParameter {
	double Mu[cNum][Dim] = { 0.0 }; // True value, Mu, Sigma are related to p(yi=0|X,theta) , p(yi=1|X,theta)
	double elnMu[cNum][Dim] = { 0.0 }; // eln value
	double Sigma[cNum][Dim][Dim] = { 1.0 };
	double Pi_orig = 0.5;

	double Epsilon = 0.01; // Epsilon: p(zi=0|zpi=1) = Epsilon, is related to p(zi=0, zpi=1|X,theta}, p(zi=1, zpi=1|X,theta) when i has parents
	double Pi = 0.5; // Pi: p(zi = 1) = Pi, is related to p(zi=0|X,theta}, p(zi =1|X,theta) when i has no parents
					 //double M = 0.4; // M: M = y0GivenZ1, is related to p(yi=0,zi=1|X,theta), p(yi=1, zi=1|X,theta)
	int CollapseSwitch = 0;
	double rho = 0.5;

	vector<double> elnPxn_zn;
	vector<double> elnPzn_xn;
	double elnPzn[cNum];
	double elnPy_z[cNum][cNum]; //transition probabilities
	double elnPz[cNum];
	double elnPz_zpn[cNum][cNum];

	int allPixelSize = -1; //set a wrong default value, so it will report error if it's not initialized properly
	double THRESHOLD = 0.05;
	int maxIteratTimes = 30;
	int maxParentDegree = 10;
	int ROW = -1;
	int COLUMN = -1;
	int NAValue = -1;
	int orgPixelSize = -1;
	int bfsRoot = -1;
	string reachId;
	string fname;
	float cutoff_percentage=0.2;
	float minCost = 0.5;
	int useCutoff = 1; // Added by Saugat: parameter to decide whether to use cutoff or not in FIST
	int useHMT = 1; // whether or not to use HMT tree (SplitTree)
	double lambda;
	
	ostringstream lambda_str;
    
};

// Added by Saugat: for multi reach inference
struct mData {
	//vector<pair<int, vector<int, float>>> maxCost_nodeid_Pairs; // reach_num(13,42,49..), 
};


struct sData {
	vector<bool>NA;
	vector<float>features; // RGB +..., rowwise, long array
	vector<double>elevationVector;
	vector<int> sortedElevationIndex;
	vector< pair<double, int> >elevationIndexPair; // here index means expanded index
	vector<Node*>allNodes;
	vector<Node*>leftallNodes;
	vector<Node*>rightallNodes;
	vector<double>A;
	vector<double>B;
	vector<int>C;
	vector<int> bfsTraversalOrder;
	vector<vector<int>> leftbfsOrder;
	vector<vector<int>> rightbfsOrder;
	vector<int> leftbfsRootNodes;
	vector<int> rightbfsRootNodes;
	vector<int> AdjustedReachNodes;
	vector<double> highestCostLeft;
	vector<double> highestCostRight;
	vector<double> inferedmaxCostLeft;
	vector<double> inferedmaxCostRight;
	vector<double> combinedCost;
	vector<double> avgCost;
	vector<int> hasObservedPixelsLeft;
	vector<int> hasObservedPixelsRight;
	vector<int> reach_ids;
	vector<int> reach_ids_orig;
	vector<int> river_ids;
	vector<int> river_ids_orig;

	map<int, bool> river_ids_map;
	map<int, bool> reach_ids_map;
	map<int, bool> reach_ids_orig_map;
	map<int, bool> reaches_map;

	vector<int> reach_ids_skipped_orig;
	vector<int> reach_ids_skipped;
	vector<int> reaches;
	vector<int> distance_covered;
	vector<double> flow_accumulation;
	vector<double> reach_fel;
	vector<int> leftbfsRootNodes_orgId;
	vector<int> rightbfsRootNodes_orgId;
	
	//for testing purpose
	vector<double> inferedcostLeft_afterInference;
	vector<double> inferedcostRight_afterInference;
	vector<double> inferedcostLeft_afterBranchStats;
	vector<double> inferedcostRight_afterBranchStats;
	vector<double> inferedcostLeft_afterInterpolation;
	vector<double> inferedcostRight_afterInterpolation;

	vector<pair<int, float>> chainLength_cost_Pairs;
	vector<pair<float, float>>maxChainCost_cost_Pairs;

	vector<pair<float, int>> maxCost_nodeid_Pairs;
	vector<pair<int, int>>chainLength_nodeid_Pairs;
	vector<int> boundary_LeafNodes;
	map<int, int> leafNode_boundaryNodes;
	std::set<pair<float, int>> cost_boundaryNode_Pairs;
	map<int, vector<int>> boundaryNode_leafNodes_Map;

	// TODO: added by Saugat
	vector<vector<double>> costVectorLeft;
	vector<vector<int>> sortedCostIndexLeft;
	vector< vector<pair<double, int>> > costIndexPairLeft; // here index means expanded index

	vector<vector<double>> costVectorRight;
	vector<vector<int>> sortedCostIndexRight;
	vector< vector<pair<double, int>> > costIndexPairRight; // here index means expanded index
	vector<int> intervalVal;
	vector < unordered_map<size_t, int> > leftIndex2PixelId;
	vector < unordered_map<size_t, int> > rightIndex2PixelId;

	//vector<double> * costVectorLeft;
	//vector<int> * sortedCostIndexLeft;
	//vector< pair<double, int> > * costIndexPairLeft; // here index means expanded index

	//vector<double>* costVectorRight;
	//vector<int>* sortedCostIndexRight;
	//vector< pair<double, int> >* costIndexPairRight; // here index means expanded index

	// added by Saugat
	vector<int> leftInferredRegionsOld;
	//vector<int> rightInferredRegions;

	map<int, bool> leftInferredRegions;
	map<int, bool> rightInferredRegions;
	map<int, double> inferredFloodFrontier_regionId2Cost;
	map<int, double> inferredFloodFrontier_regionId2nodeId; 
	map<int, size_t> inferredFloodFrontier_regionId2Size;
// 	vector<double> inferredFloodFrontier;
	
	map<int, bool> leftBankNodes;
	
	vector<int> leftNodesInOrder;
	vector<int> rightNodesInOrder;
	vector<int> leftNodesInCorrectOrder;
	vector<int> rightNodesInCorrectOrder;
	
	vector<double> inferedmaxCostRight_new;
	vector<double> inferedmaxCostLeft_new;
	
	map<int, double> inferedmaxCost_id2cost;
// 	map<int, double> inferedmaxCostRight_id2cost;

    map<int, int> inferedmaxCostLeft_idx2id;
    map<int, int> inferedmaxCostRight_idx2id;


	vector<bool> specialDataLeft;
	vector<bool> specialDataRight;
	// for loglikelihood
	vector<vector<double>> loglikelihood_leftRegions;
	vector<vector<double>> loglikelihood_rightRegions;

	vector<vector<int>> mainBranchNodeIds_leftRegions;
	vector<vector<int>> mainBranchNodeIds_rightRegions;

	// for viterbi
	vector<size_t> result_ids_left;
	vector<double> result_costs;

	vector<size_t> result_ids_right;
	vector<double> result_costs_right;
	
	map<int, bool> large_regions;
	
	vector<double> regularizedMaxCostLeft;
	vector<double> regularizedMaxCostRight;
	
	map<int, int> leftOrderSelectedToLeftOrder; // map for selected left order to actual left order
	map<int, int> rightOrderSelectedToRightOrder; // map for selected left order to actual left order



};


struct sTree {
	queue<int> leavesIDQueue;
	vector<int>nodeLevels;
};


struct sInference {
	vector<double> incoming;
	vector<double> outgoing;

	vector<double>waterLikelihood;  //-0.5(x-mu)T Sigma^-1 (x-mu) for class 1
	vector<double>dryLikelihood;  //-0.5(x-mu)T Sigma^-1 (x-mu) for class 0
	double lnCoefficient[cNum]; //log constanct of Gaussian distribution for two classes

	vector<double> lnTop2y;
	vector<double> lnz2top; //from z to y, outgoing of z, before factor node
	vector<double> lnbottom2y; //Vertical incoming message from z to y, after factor node

	vector<double> lnvi; // Top2z message
	vector<double> lnfi; // leaves to root, incoming messages; from zp to z except zp=NULL, after factor node
	vector<double> lnfo; // leaves to root, outgoing messages; product of lnfi and lngi, before factor node
	vector<double> lngi; // root to leaves, incoming messages; Horizontal incoming message from z to zp, after factor node
	vector<double> lngo; // root to leaves, outgoing messages

	vector<double> marginal_Yn;
	vector<double> marginal_YnZn; // only zn = 1 is useful to update parameters
	vector<double> marginal_ZnZpn; // only Zpn = 1 is useful to update parameters
	vector<double> marginal_Zn; // Separate from ZnZpn now. Those internal values are just for testing.

	vector<int>chainMaxGainFrontierIndex;

	//added by Arpan
	vector<vector<int>> correspondingNeighbours;
	vector<vector<int>> correspondingNeighboursClass;
};

struct subset
{
	int parent;
	int rank;
};


struct conMatrix {
	int TT;
	int TF;
	int FF;
	int FT;
};


// Added by Saugat
struct RegionNode {
	double regionCenterRow;
	double regionCenterCol;
	double regionLoc; // how far is region from reach center
	bool regionTooFar;

	RegionNode() {

	}
};

struct extras {
	double reachLength;
	double reachCenterRow;
	double reachCenterCol;

	vector<RegionNode*>leftRegions;
	vector<RegionNode*>rightRegions;

	vector<double> standardDeviationLeft;
	vector<double> standardDeviationRight;
};

// struct sRegularization{
// 	queue<pair<int, int>> left_bfs_que;
// 	map<int, bool> left_bfs_visited;
// 	vector<int> left_node_order;
// 	vector<int> right_node_order;
// 	map<int, bool> on_queue_left;
// };




