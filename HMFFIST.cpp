//Assumption no missing Value in this test data
#define Dim 1 //input data dimension
#define cNum 2
#define _USE_MATH_DEFINES
#define LOGZERO -INFINITY
#define LOGLARGE -1000000
//#define LOGZERO -INFINITY
#define _CRT_SECURE_NO_WARNINGS
#define MESSAGELOW -INFINITY
#define MESSAGEHIGH 10000
#define PIXELLIMT 25000
#define MAXCOST -1000.0
#define MAXGAIN INFINITY

#define EB 2  //1 use max cost for identifying effective branches 0 use chain length
#define BOUNDARY_NODES_OBSERVED 2   // 0 does not exclude any branches // 1 consider pits and tree and unobserved and exclude those branches // 2 excludes branches if the boundary of flood and dry is overlapping with pits layer
#define EFFECTIVE_BRANCH_VIZ 1
#define DEBUG_OUTPUT 1
#include<iostream>
#include<functional>
#include <sys/types.h>
#include <sys/stat.h>
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
#include<unordered_set>
#include <iomanip>
#include <sstream>
#include <map>
#include "GeotiffRead.cpp"
#include "GeotiffWrite.cpp"
//#include "Tree.cpp"
#include "DataTypes.h"


using namespace std;
double Pi_orig = 0.5;
ofstream timeSave;
float SUM_COST = 0.0;
int COUNTER = 0;

std::map<int, int> nodeIndexMap;

bool comp(Node* a, Node* b);
bool comp(Node* a, Node* b)
{
	return a->cost < b->cost;
}
class cFlood {
private:
	//struct sComponent unionComponent;
	struct sParameter parameter;
	struct sData data;
	vector<struct sData> subtreeData;
	struct sTree tree;
	struct sInference infer;
	vector<int>testIndex;
	vector<int>testLabel;
	vector<int>mappredictions;
	//ofstream timeLogger;
	std::string CTInputLocation;
	std::string CTSourceDirection;
	std::string CTProbability;
	std::string CTBank;
	std::string CTCost;
	std::string CTPits;
	std::string CTTree;
	std::string CTRoughness;
	std::string CTFel;
	std::string CTPara;
	std::string CTStream;
	std::string(CTOutputFolderByDate);
	std::string CTOutputLocation;
	std::string CTPrediction;
	std::string CTLeftBank;
	std::string CTPredictionTxt;
	std::string CTPredictionRegularized;
	std::string CTPredictionRegularizedTxt;
	std::string CTTestCase;
	std::string CTAlgo;

	//tree construction
	struct subset* subsets;

	//new
	vector<double>elnPzn_xn;

	// added by Saugat
	struct extras extra;

public:
	void input(int argc, char* argv[]);


	void UpdateTransProb(); //Update P(y|z), P(zn|zpn), P(zn|zpn=empty)
	void UpdatePX_Z();

	// learning
	void learning();
	void MessagePropagation();
	void UpdateMarginalProb();
	void UpdateParameters();


	//inference
	void inference();
	void interpolate();
	void output();
	void prediction();
	void prediction_FIST();
	void selected_prediction();
	void selected_prediction_FIST();

	//helper functions
	void removeLink(vector<int>& v, int removeID);
	void displayTree(int TreeID);
	void updateMapPrediction_left();
	void updateMapPrediction_right();

	vector<int> getBFSOrder(int root, vector<int>& bfsVisited, int bank);
	//struct conMatrix getConfusionMatrix();

	void updateMapPrediction_left_new();
	void updateMapPrediction_left_hmt();
	void updateMapPrediction_right_hmt();
	void updateMapPrediction_right_verify();

	void verify_deltaResult_left();
	void verify_deltaResult_right();
	void delta_prediction();

	//Test Modules
	void sanityChecker();
	void getOriginIdBanks();
	void getOrgIds();
	void getIds();
	void getOriginIdLeftBanks();
	void getOriginIdRightBanks();
	void getRegionNodeCount();
	void getLeftRegionNodeCount();
	void getRightRegionNodeCount();
	void getOriginIdBanks_effectiveBranches();

	// -- added by Saugat --
	// hmt tree
	void splitTree();
	void getNewBFSOrder();

	//utilities
	int find(struct subset subsets[], int i);
	void Union(struct subset subsets[], int x, int y);
	void validateTreeLeft();
	void validateTreeRight();

	void validateTreeInferenceLeft();
	void validateTreeInferenceRight();

	void validateTreeInferenceLeftFIST();
	void validateTreeInferenceRightFIST();

	void getStatistics();

	void export_FIST_structure(string child_file, string parent_file, string small_file);

	// get left and right nodes order
// 	void getNodeOrder(queue<pair<int, int>>& bfs_que, map<int, bool>& bfs_visited, vector<int>& left_node_order, map<int, bool>& on_queue);
	int reachBFS(queue<pair<int, int>>& bfs_que, map<int, bool>& bfs_visited, vector<int>& left_node_order, map<int, bool>& on_queue); // for black nodes
// 	void brokenBFS(int last_node, queue<pair<int, int>>& bfs_que, map<int, bool>& bfs_visited, map<int, bool>& on_queue); // after the chain is broken

	// log likelihood regularization
	void getLoglikelihood();

	// viterbi algorithm
	void viterbi(double lambda, int ranges);

	void saveExtraInfo();
	// for later start------
	void selected_prediction_regularized(string lambda);
	void prediction_regularized();
	void interpolate_regularized(string lambda);
	void outputRegulization(string lambda);
	// clear all the vectors
	void clear_all();
};



void getCofactor(double mat[Dim][Dim], double temp[Dim][Dim], int p, int q, int n) {
	int i = 0, j = 0;
	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				temp[i][j++] = mat[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

//dynamic memory allocation,dimensional two dimension array
/* Recursive function for finding determinant of matrix.
n is current dimension of mat[][]. */
double determinant(double mat[Dim][Dim], int n) {
	double D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1)
		return mat[0][0];

	double temp[Dim][Dim]; // To store cofactors
	int sign = 1;  // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++) {
		// Getting Cofactor of mat[0][f]
		getCofactor(mat, temp, 0, f, n);
		D += sign * mat[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}
	return D;
}

void adjoint(double A[Dim][Dim], double adj[Dim][Dim]) {
	if (Dim == 1) {
		adj[0][0] = 1;
		return;
	}

	int sign = 1;
	double temp[Dim][Dim];

	for (int i = 0; i < Dim; i++) {
		for (int j = 0; j < Dim; j++) {
			getCofactor(A, temp, i, j, Dim);


			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(temp, Dim - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(double A[Dim][Dim], double inverse[Dim][Dim]) {
	// Find determinant of A[][]

	if (Dim == 1) {
		inverse[0][0] = 1.0 / A[0][0];
		return true;
	}

	double det = determinant(A, Dim);
	if (det == 0) {
		std::cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	double adj[Dim][Dim];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			inverse[i][j] = adj[i][j] / double(det);
	return true;
}

// extended ln functions
double eexp(double x) {
	if (x == LOGZERO) {
		return 0;
	}
	else {
		return exp(x);
	}
}

double eln(double x) {
	if (x == 0) {
		return LOGZERO;
	}
	else if (x > 0) {
		return log(x);
	}
	else {
		std::cout << "Negative input error " << x << endl;
		exit(0);
	}
}

double eln_ll(double x) {
	if (x == 0) {
		return LOGLARGE;
	}
	else if (x > 0) {
		return log(x);
	}
	else {
		std::cout << "Negative input error " << x << endl;
		exit(0);
	}
}

double elnsum(double x, double y) {
	if (x == LOGZERO) {
		return y;
	}
	else if (y == LOGZERO) {
		return x;
	}
	else if (x > y) {
		return x + eln(1 + eexp(y - x));
	}
	else {
		return y + eln(1 + eexp(x - y));
	}
}

double elnproduct(double x, double y) {
	if (x == LOGZERO || y == LOGZERO) {
		return LOGZERO;
	}
	else {
		return x + y;
	}
}
int dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
	// return (info.st_mode) ? 1 : 0;
}

bool dirExists_2(const char* const s) {
	struct stat buffer;
	return (stat(s, &buffer) == 0);
}


void cFlood::UpdateTransProb() {
	if (cNum != 2) {
		std::cout << "cannot handle more than two classes now!" << endl;
		std::exit(1);
	}

	double eln(double);
	parameter.elnPz[0] = eln(1 - eexp(parameter.Pi));
	parameter.elnPz[1] = parameter.Pi;
	parameter.elnPz_zpn[0][0] = eln(1);
	parameter.elnPz_zpn[0][1] = parameter.Epsilon;
	parameter.elnPz_zpn[1][0] = eln(0);
	parameter.elnPz_zpn[1][1] = eln(1 - eexp(parameter.Epsilon));
	if (eexp(parameter.Epsilon) < 0 || eexp(parameter.Epsilon) > 1) {
		std::cout << "Epsilon Error: " << eexp(parameter.Epsilon) << endl;
	}
	if (eexp(parameter.Pi) < 0 || eexp(parameter.Pi) > 1) {
		std::cout << "Pi Error: " << eexp(parameter.Pi) << endl;
	}
	if (eexp(parameter.elnPz_zpn[0][1]) + eexp(parameter.elnPz_zpn[1][1]) != 1) {
		std::cout << "Error computing parameter.elnPz_zpn " << endl;
	}
	if (eexp(parameter.elnPz[0]) + eexp(parameter.elnPz[1]) != 1) {
		std::cout << "Error computing parameter.elnPz " << endl;
	}
}

void cFlood::UpdatePX_Z() {
	// Calculate inverse of sigma
	double adjointMatrix[cNum][Dim][Dim]; // To store adjoint of A[][]
	double inverseMatrix[cNum][Dim][Dim]; // To store inverse of A[][]
	for (int c = 0; c < cNum; c++) {
		adjoint(parameter.Sigma[c], adjointMatrix[c]);
	}
	for (int c = 0; c < cNum; c++) {
		if (!inverse(parameter.Sigma[c], inverseMatrix[c])) {
			cout << "Inverse error" << endl;
		}
	}

	//xiGivenZi_coefficient, log form
	for (int c = 0; c < cNum; c++) {// |Sigma|^(-1/2)
		infer.lnCoefficient[c] = -0.5 * Dim * log(2 * M_PI) - 0.5 * log(fabs(determinant(parameter.Sigma[c], Dim)));
	}

	// Calculate p(x|z)
	double intermediateValue[cNum][Dim] = { 0 };
	double likelihood[cNum] = { 0 };
	double xMinusMu[cNum][Dim] = { 0 };

	for (size_t i = 0; i < parameter.allPixelSize; i++) {
		if (!data.NA[i]) { // Not missing data

			for (int c = 0; c < cNum; c++) {
				likelihood[c] = 0;
			}

			for (int c = 0; c < cNum; c++) {
				for (int d = 0; d < Dim; d++) {
					intermediateValue[c][d] = 0;
				}
			}

			// -0.5*(x-mu)' * Sigma^-1 * (x-mu), matrix multiply
			for (int c = 0; c < cNum; c++) {
				for (int d = 0; d < Dim; d++) {
					xMinusMu[c][d] = data.features[i * Dim + d] - parameter.Mu[c][d];
				}
			}

			for (int c = 0; c < cNum; c++) {
				for (int k = 0; k < Dim; k++) {
					for (int n = 0; n < Dim; n++) {
						intermediateValue[c][k] += xMinusMu[c][n] * inverseMatrix[c][n][k];
					}
					likelihood[c] += intermediateValue[c][k] * xMinusMu[c][k];
				}
			}

			for (int cls = 0; cls < cNum; cls++) {
				parameter.elnPxn_zn[i * cNum + cls] = -0.5 * likelihood[cls] + infer.lnCoefficient[cls];
			}

		}
		else {
			for (int cls = 0; cls < cNum; cls++) {
				parameter.elnPxn_zn[i * cNum + cls] = eln(1);
			}
		}

	}
}
// selected prediction after regularization
void cFlood::selected_prediction_regularized(string lambda) {
	cout << "Selected prediction regularized started!" << endl;

	std::fill(mappredictions.begin(), mappredictions.end(), -1);



	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		/*if (data.hasObservedPixelsRight[rightOrder] && data.rightbfsOrder[rightOrder].size()>= PIXELLIMT) {
			continue;
		}*/
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];

				// commenting temporarily
				if (data.regularizedMaxCostRight[rightOrder] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					continue;
					/*if (data.allNodes[nodid]->isNa == 0)
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder] * 2;*/
				}
				else {

					// Refill
					if (data.allNodes[nodid]->cost <= data.regularizedMaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) { // CHECK NC
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder];
					}
					else {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}




				}
			}
		}
	}

	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		mappredictions[data.allNodes[data.rightNodesInOrder[i]]->originalId] = 1;
	}

	// nodes in river(high prob pixels from UNet should also be flooded)
	for (int i = 0; i < data.river_ids.size(); i++) {
		mappredictions[data.allNodes[data.river_ids[i]]->originalId] = 1;
	}

	// Comment: END

	auto start = std::chrono::system_clock::now();

	ofstream classout;
	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}

	string tree_type = "fist_tree";
	if (parameter.useHMT == 1) {
		tree_type = "hmt_tree";
	}

	idf = idf + "_" + tree_type;

	classout.open(CTOutputLocation + "selected_" + "ranges_" + lambda + "_" + CTPredictionRegularizedTxt);

	for (int i = 0; i < mappredictions.size(); i++) {
		classout << mappredictions[i] << endl;
	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
			prediction[row][col] = mappredictions[index];
			index++;
		}
	}
	GDALDataset* srcDataset = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform[6];
	srcDataset->GetGeoTransform(geotransform);
	const OGRSpatialReference* poSRS = srcDataset->GetSpatialRef();
	;
	GeotiffWrite finalTiff((CTOutputLocation + "selected_" + "ranges_" + lambda + "_" + CTPredictionRegularized).c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform, poSRS);
	finalTiff.write(prediction);

	cout << "Selected Prediction Regularized finished!" << endl;
}


void cFlood::prediction_regularized() {
	cout << "prediction regularized started!" << endl;


	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

		// no need to touch inferred regions
		if (data.rightInferredRegions[data.rightNodesInOrder[rightOrder]]) {
			continue;
		}

		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				if (data.regularizedMaxCostRight[rightOrder] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					continue;
				}
				else {

					if (data.allNodes[nodid]->cost <= data.regularizedMaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) { // CHECK NC
						// mappredictions[data.allNodes[nodid]->originalId] = 1;
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder];
					}
					/*else if(data.allNodes[nodid]->isNa == 0){*/
					else if (mappredictions[data.allNodes[nodid]->originalId] == 0) {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}
					else {
						mappredictions[data.allNodes[nodid]->originalId] = -1; // NC
					}
				}
			}
		}
	}

	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		mappredictions[data.allNodes[data.rightNodesInOrder[i]]->originalId] = 1;
	}


	cout << "prediction regularized finished!" << endl;
}

void cFlood::interpolate_regularized(string lambda) {
	cout << "interpolation regularized started!" << endl;
	//profile table before interpolation
	ofstream profiletable;

		//find the regions with infered cost values (i.e values that are not -1)
	vector<int> stops;

	ofstream profiletable_right;


	profiletable_right.open(CTOutputLocation + "ProfileTables/" + CTTestCase + "_ProfileTable_preInterpolation_right_ranges_" + lambda + "_Viterbi" + ".csv");
	//profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.fname + "_ProfileTable_preInterpolation.csv");
	profiletable_right << "SourceId" << "," << "Right Regularized Cost" << endl;
	for (int index = 0; index < data.rightNodesInOrder.size(); index++) {
		profiletable_right << data.rightNodesInOrder[index] << ","
			<< data.regularizedMaxCostRight[index] << endl;
	}
	profiletable_right.close();

	//for right bank.
	int current = 0;
	while (current < data.rightNodesInOrder.size()) {
		if (data.regularizedMaxCostRight[current] == -1 && current == 0) {

			//find the first reach node with non -1 max cost value
			int index = -1;
			for (int j = 1; j < data.rightNodesInOrder.size(); j++) {
				if (data.regularizedMaxCostRight[j] != -1) {
					index = j;
					break;
				}
			}
			if (index == -1) {
				break;
			}
			double value = data.regularizedMaxCostRight[index];
			for (int i = 0; i < index; i++) {
				data.regularizedMaxCostRight[i] = value;
				data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = value;
			}
			current = index;


		}
		else if (data.regularizedMaxCostRight[current] != -1) {
			//two cases
				//case 1: there are n points in between next reach that has cost value
				//case 2: there is no next point
			//find index of next reach node that has cost value
			int index = -1;
			int count = 0;
			double value = data.regularizedMaxCostRight[current];
			for (int j = current + 1; j < data.rightNodesInOrder.size(); j++) {
				if (data.regularizedMaxCostRight[j] != -1) {
					index = j;
					break;
				}
				count++;
			}
			if (index == -1) {// case 2
				for (int i = current + 1; i < data.rightNodesInOrder.size(); i++) {
					data.regularizedMaxCostRight[i] = value;
					data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = value;
				}
				current = data.rightNodesInOrder.size();
				break;
			}
			else if (count == 0 && index == current + 1) {
				current = index;
			}
			else {
				double interval = (data.regularizedMaxCostRight[index] - value) / count;
				for (int i = current + 1; i < index; i++) {
					data.regularizedMaxCostRight[i] = data.regularizedMaxCostRight[(i - 1)] + interval;
					data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = data.regularizedMaxCostRight[(i - 1)] + interval;
				}
				current = index;
			}
		}

	}
	cout << "interpolation regularized finished!" << endl;

}
void cFlood::outputRegulization(string lambda) {
	auto start = std::chrono::system_clock::now();

	ofstream classout;
	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}

	string tree_type = "fist_tree";
	if (parameter.useHMT == 1) {
		tree_type = "hmt_tree";
	}

	idf = idf + "_" + tree_type;

	classout.open(CTOutputLocation + CTPredictionRegularizedTxt);


	for (int i = 0; i < mappredictions.size(); i++) {
		int val = 0;
		if (mappredictions[i] > 0) {
			val = 1;
		}
		classout << val << endl;

	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
		    if(mappredictions[index] >0)
		        prediction[row][col] = 1;
		    else
			    prediction[row][col] = mappredictions[index];
			index++;

		}
	}
	GDALDataset* srcDataset = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform[6];
	srcDataset->GetGeoTransform(geotransform);
	const OGRSpatialReference* poSRS = srcDataset->GetSpatialRef();

	GeotiffWrite finalTiff((CTOutputLocation + CTPredictionRegularized).c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform, poSRS);
	finalTiff.write(prediction);



	ofstream profiletable;
	profiletable.open(CTOutputLocation + "ProfileTables/" + CTTestCase + "_ProfileTable_right_ranges_" + lambda + "_Viterbi" + ".csv");
	//profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.fname + "_ProfileTable.csv");
	profiletable << "SourceId" << "," << "Right Regularized Cost" << endl;
	for (int index = 0; index < data.rightNodesInOrder.size(); index++) {
		profiletable
			<< data.rightNodesInOrder[index] << ","
			<< data.inferedmaxCostRight[index] << endl;
	}
	profiletable.close();

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double>elapsed_seconds = end - start;
	std::cout << "Writing Prediction File took " << elapsed_seconds.count() << "seconds" << endl;

	//clear_all();
}
//Assume the first node is node without parents
//Assume the first node is node without parents
void cFlood::MessagePropagation() {
	//NOTE: we can only handle 64 parents/children for long int bit_counter;
	//sort all nodes in BFS traversal order
	vector<int> mpVisited(parameter.allPixelSize, 0);
	//leaves to root

	// bfsTraversalOrder

	// go through each regions
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) { //// go through every reach ids

		int bfsTraversalOrderSize = (int)data.rightbfsOrder[rightOrder].size();

		for (int node = bfsTraversalOrderSize - 1; node >= 0; node--) {
			int cur_node_id = data.rightbfsOrder[rightOrder][node];  //n
			//initializing fi_childlist,fi fo
			data.allNodes[cur_node_id]->fi_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[cur_node_id]->fi[cls] = 0;
				data.allNodes[cur_node_id]->fo[cls] = 0;
			}

			//first figure out which neighbor fmessage passes to from current node pass n->? foNode;
			//idea: In bfs traversal order leave to root, check if next the node in bfs order is parent or child of the current node (should be child or parent of the current node)
			int foNode = -1;
			bool foNode_isChild = false;
			for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {  //check in parent list if found respective parent node is foNode
				int pid = data.allNodes[cur_node_id]->parentsID[p];
				if (!mpVisited[pid]) {
					foNode = pid;
					break;
				}
			}
			if (foNode == -1) {
				for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
					int cid = data.allNodes[cur_node_id]->childrenID[c];
					if (!mpVisited[cid]) {
						foNode = cid;
						foNode_isChild = true;
						break;
					}
				}
			}
			data.allNodes[cur_node_id]->foNode = foNode;
			data.allNodes[cur_node_id]->foNode_ischild = foNode_isChild;
			int root_nodeId = data.rightbfsRootNodes[rightOrder];
			if (cur_node_id == root_nodeId && data.allNodes[cur_node_id]->childrenID.size() == 0) {
				foNode_isChild = true;
				data.allNodes[cur_node_id]->foNode_ischild = true;
			}
			//incoming message from visited child
			if (data.allNodes[cur_node_id]->childrenID.size() > 0) {

				for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
					int child_id = data.allNodes[cur_node_id]->childrenID[c];

					if (child_id == foNode) {  //if child_node is foNode skip

						continue;
					}
					//extract parents except current node
					vector<int> parentOfChildExceptCurrentNode;   //Yk E Pc k!=n
					for (int en = 0; en < data.allNodes[child_id]->parentsID.size(); en++) {
						if (data.allNodes[child_id]->parentsID[en] == cur_node_id) { //k!=n
							continue;
						}
						parentOfChildExceptCurrentNode.push_back(data.allNodes[child_id]->parentsID[en]);
					}

					for (int cls = 0; cls < cNum; cls++) {  //cls represents current node class
						double sumAccumulator = eln(0);   //should be 0 since we are summing it up//eln(1);//need to confirm
						for (int c_cls = 0; c_cls < cNum; c_cls++) { //c_cls reperesnets child class label   Yc
							int max_bitCount = 1 << parentOfChildExceptCurrentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent and child class label(given by c_cls)
								double productAccumulator = eln(1);
								int parentClsProd = 1; //p(c), product of parent classes for child c

								for (int p = 0; p < parentOfChildExceptCurrentNode.size(); p++) {//calculating Product(fo(p)) for all parent of current child except the current node
									int pid = parentOfChildExceptCurrentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
								}
								productAccumulator = elnproduct(productAccumulator, data.allNodes[child_id]->fo[c_cls]);  //product with fo(c)
								//multiplying P(Yc|Ypc)
								parentClsProd *= cls; //class of current node
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[c_cls][parentClsProd]);
								cout << "Line 508" << endl;
								sumAccumulator = elnsum(sumAccumulator, productAccumulator);
							}
						}
						data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls] = sumAccumulator;
					}
				}
			}

			if (foNode_isChild) {  //means the current node has all visited parents
				if (data.allNodes[cur_node_id]->parentsID.size() == 0) {
					for (int cls = 0; cls < cNum; cls++) {
						data.allNodes[cur_node_id]->fi[cls] = parameter.elnPz[cls];
					}
				}
				else {
					for (int cls = 0; cls < cNum; cls++) {
						double sumAccumulator = eln(0);
						int max_bitCount = 1 << data.allNodes[cur_node_id]->parentsID.size();
						for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label
							double productAccumulator = eln(1);
							int parentClsProd = 1;
							for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
								int pid = data.allNodes[cur_node_id]->parentsID[p];
								int parentClsValue = (bitCount >> p) & 1;
								parentClsProd *= parentClsValue;
								productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
							}
							productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[cls][parentClsProd]);
							cout << "Line 537" << endl;
							sumAccumulator = elnsum(sumAccumulator, productAccumulator);
						}
						data.allNodes[cur_node_id]->fi[cls] = sumAccumulator;
					}
				}

				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
						int child_id = data.allNodes[cur_node_id]->childrenID[c];

						if (child_id == foNode) {
							continue;
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls]); //multiplying with al the child fi except the outgoing child
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi[cls]);  // multiplying with fi(n)_parent
					productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}

			}

			else {  //message pass n-> parent there is no fi(n)_parent   //computes for root node as well
				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls]); //multiplying with al the child fi except the outgoing child
					}
					productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}
			}

			mpVisited[cur_node_id] = 1;

			//verification
			for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				if (data.allNodes[cur_node_id]->childrenID[c] == foNode) {
					for (int cls = 0; cls < cNum; cls++) {
						if (data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls] != 0) {
							cout << " fi_childlist Message Computation Error (this should not be computed)  Node " << cur_node_id << endl;
						}
					}
					continue;
				}
				for (int cls = 0; cls < cNum; cls++) {
					if (data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls] < MESSAGELOW || data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls]>0) {
						cout << " fi_childlist Message Computation Error in Node " << cur_node_id << endl;
					}
				}
			}
			if (foNode_isChild) {
				for (int cls = 0; cls < cNum; cls++) {
					if (data.allNodes[cur_node_id]->fi[cls] < MESSAGELOW || data.allNodes[cur_node_id]->fi[cls]>0) {
						cout << " fi Message Computation Error in Node " << cur_node_id << endl;
					}
				}
			}
			//verify fo
			for (int cls = 0; cls < cNum; cls++) {
				if (data.allNodes[cur_node_id]->fo[cls] < MESSAGELOW || data.allNodes[cur_node_id]->fo[cls]>0) {
					cout << " fo Message Computation Error in Node " << cur_node_id << endl;
				}
			}

		}


		//root to leaves traversal
		vector<int> gVisited(parameter.allPixelSize, 0);
		//for root node
		//computing gi
		//int root_nodeId = data.bfsTraversalOrder[0]; //root node
		//cout << "root_nodeId = " << root_nodeId << endl;

		int root_nodeId = data.rightbfsRootNodes[rightOrder];
		if (data.allNodes[root_nodeId]->childrenID.size() == 0) {
			//if (data.allNodes[root_nodeId]->childrenID.size() == 0) {
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[root_nodeId]->go_parent[cls] = 0;
				data.allNodes[root_nodeId]->gi[cls] = eln(1);
			}
			data.allNodes[root_nodeId]->go_fromParent = false;        //case1: if current node has visited parent
			data.allNodes[root_nodeId]->go_fromChild = true;         //case2: if current node has visited child and if go is in visited child
			data.allNodes[root_nodeId]->go_fromParentofChild = false; //case3: if current node has visited child and if go is in one of the visited child's parent
			data.allNodes[root_nodeId]->go_lastVisitedNode = -1;
			for (int cls = 0; cls < cNum; cls++) {
				double productAccumulator = eln(1);
				productAccumulator = elnproduct(productAccumulator, data.allNodes[root_nodeId]->gi[cls]);
				productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[root_nodeId * cNum + cls]);
				data.allNodes[root_nodeId]->go_parent[cls] = productAccumulator;
			}

		}
		else {
			//initializing go_childlist,go_parent gi for root node
			data.allNodes[root_nodeId]->go_ChildList.resize(data.allNodes[root_nodeId]->childrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[root_nodeId]->go_parent[cls] = 0;
				data.allNodes[root_nodeId]->gi[cls] = parameter.elnPz[cls];
			}
			data.allNodes[root_nodeId]->go_fromParent = true;        //case1: if current node has visited parent
			data.allNodes[root_nodeId]->go_fromChild = false;         //case2: if current node has visited child and if go is in visited child
			data.allNodes[root_nodeId]->go_fromParentofChild = false; //case3: if current node has visited child and if go is in one of the visited child's parent
			data.allNodes[root_nodeId]->go_lastVisitedNode = -1;
			//computing go for every child c of n
			for (int c = 0; c < data.allNodes[root_nodeId]->childrenID.size(); c++) {
				int cid = data.allNodes[root_nodeId]->childrenID[c];
				for (int cls = 0; cls < cNum; cls++) {
					double productAccumulator = eln(1);
					for (int d = 0; d < data.allNodes[root_nodeId]->childrenID.size(); d++) {
						if (d == c) continue;
						productAccumulator = elnproduct(productAccumulator, data.allNodes[root_nodeId]->fi_ChildList[d * cNum + cls]);
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[root_nodeId]->gi[cls]);
					productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[root_nodeId * cNum + cls]);
					data.allNodes[root_nodeId]->go_ChildList[c * cNum + cls] = productAccumulator;
				}
			}
		}
		gVisited[root_nodeId] = 1;
		for (int node = 1; node < data.rightbfsOrder[rightOrder].size(); node++) {
			//for (int node = 1; node < data.bfsTraversalOrder.size(); node++) {
			int cur_node_id = data.rightbfsOrder[rightOrder][node];
			//int cur_node_id = data.bfsTraversalOrder[node];  //n
															 //only one gi, many go
															 //initializing go_childlist,go_parent gi
			data.allNodes[cur_node_id]->go_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[cur_node_id]->go_parent[cls] = 0;
				data.allNodes[cur_node_id]->gi[cls] = 0;
			}
			//data.allNodes[cur_node_id]->go_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size()*cNum, LOGZERO);

			//first figure out g direction from parent side or      child side (two case: child side or parent of child side)
			data.allNodes[cur_node_id]->go_fromParent = false;        //case1: if current node has visited parent
			data.allNodes[cur_node_id]->go_fromChild = false;         //case2: if current node has visited child and if go is in visited child
			data.allNodes[cur_node_id]->go_fromParentofChild = false; //case3: if current node has visited child and if go is in one of the visited child's parent
			data.allNodes[cur_node_id]->go_lastVisitedNode = -1;



			int visitedCounter = 0;
			for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {  //check in parent list if found respective parent node is foNode
				int pid = data.allNodes[cur_node_id]->parentsID[p];
				if (gVisited[pid]) {
					data.allNodes[cur_node_id]->go_fromParent = true;
					data.allNodes[cur_node_id]->go_lastVisitedNode = pid;
					visitedCounter++;
					break;
				}
			}
			if (data.allNodes[cur_node_id]->go_lastVisitedNode == -1) {
				for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
					int cid = data.allNodes[cur_node_id]->childrenID[c];
					if (gVisited[cid]) {
						visitedCounter++;
						if (data.allNodes[cid]->go_fromParent) {
							data.allNodes[cur_node_id]->go_fromParentofChild = true;
						}
						else {
							data.allNodes[cur_node_id]->go_fromChild = true;
						}
						data.allNodes[cur_node_id]->go_lastVisitedNode = cid;
						break;
					}
				}
			}
			if (visitedCounter != 1) {
				cout << "Not one visited Neighbour Error" << endl;
			}
			if (data.allNodes[cur_node_id]->go_fromParent == false && data.allNodes[cur_node_id]->go_fromParentofChild == false && data.allNodes[cur_node_id]->go_fromChild == false) {
				cout << "Error all neighbours are not visited Node " << cur_node_id << endl;
			}



			if (data.allNodes[cur_node_id]->go_fromParent) {

				//computing gi
				for (int cls = 0; cls < cNum; cls++) {
					double sumAccumulator = eln(0);
					int max_bitCount = 1 << data.allNodes[cur_node_id]->parentsID.size();
					for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label
						double productAccumulator = eln(1);
						int parentClsProd = 1;
						for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
							int pid = data.allNodes[cur_node_id]->parentsID[p];
							int parentClsValue = (bitCount >> p) & 1;
							parentClsProd *= parentClsValue;
							//multiply with go(Po)_childlist[n]
							if (pid == data.allNodes[cur_node_id]->go_lastVisitedNode) {
								for (int c = 0; c < data.allNodes[pid]->childrenID.size(); c++) {
									int cid = data.allNodes[pid]->childrenID[c];
									if (cid == cur_node_id) {
										double tempgoChild = data.allNodes[pid]->go_ChildList[c * cNum + parentClsValue];
										productAccumulator = elnproduct(productAccumulator, tempgoChild);
										break;
									}
								}
							}
							else {
								productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
							}
						}
						productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[cls][parentClsProd]);
						cout << "Line 747" << endl;
						sumAccumulator = elnsum(sumAccumulator, productAccumulator);
					}
					data.allNodes[cur_node_id]->gi[cls] = sumAccumulator;
				}

				//computing go for every child c of n
				for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
					int cid = data.allNodes[cur_node_id]->childrenID[c];
					for (int cls = 0; cls < cNum; cls++) {
						double productAccumulator = eln(1);
						for (int d = 0; d < data.allNodes[cur_node_id]->childrenID.size(); d++) {
							if (d == c) continue;
							productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[d * cNum + cls]);
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->gi[cls]);
						productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
						data.allNodes[cur_node_id]->go_ChildList[c * cNum + cls] = productAccumulator;
					}
				}
			}
			else {

				if (data.allNodes[cur_node_id]->go_fromChild) {
					//computing gi(n)
					int Co = data.allNodes[cur_node_id]->go_lastVisitedNode;
					vector<int> parentOfCoExcept_currentNode;
					for (int en = 0; en < data.allNodes[Co]->parentsID.size(); en++) {
						if (data.allNodes[Co]->parentsID[en] == cur_node_id) {
							continue;
						}
						parentOfCoExcept_currentNode.push_back(data.allNodes[Co]->parentsID[en]);
					}
					for (int cls = 0; cls < cNum; cls++) {  //current node class
						double sumAccumulator = eln(0);
						for (int Co_cls = 0; Co_cls < cNum; Co_cls++) {
							int max_bitCount = 1 << parentOfCoExcept_currentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label product(fo(p)) except current node
								double productAccumulator = data.allNodes[Co]->go_parent[Co_cls];
								int parentClsProd = 1;
								for (int p = 0; p < parentOfCoExcept_currentNode.size(); p++) {
									int pid = parentOfCoExcept_currentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
								}
								//p(Yco|Ypco)
								parentClsProd *= cls;
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[Co_cls][parentClsProd]);
								cout << "Line 796" << endl;
								sumAccumulator = elnsum(sumAccumulator, productAccumulator);
							}
						}
						data.allNodes[cur_node_id]->gi[cls] = sumAccumulator;
					}
				}


				else if (data.allNodes[cur_node_id]->go_fromParentofChild) {
					//computing gi(n)
					int Co = data.allNodes[cur_node_id]->go_lastVisitedNode;
					int Po = data.allNodes[Co]->go_lastVisitedNode;
					if (Po == -1) {
						cout << "error: Three should be a parent of a child" << endl;
					}
					int CIndex = -1;
					for (int c = 0; c < data.allNodes[Po]->childrenID.size(); c++) {
						if (data.allNodes[Po]->childrenID[c] == Co) {
							CIndex = c;
							break;
						}
					}
					vector<int> parentOfCoExcept_currentNode;
					for (int en = 0; en < data.allNodes[Co]->parentsID.size(); en++) {
						if (data.allNodes[Co]->parentsID[en] == cur_node_id) {
							continue;
						}
						parentOfCoExcept_currentNode.push_back(data.allNodes[Co]->parentsID[en]);
					}
					for (int cls = 0; cls < cNum; cls++) {  //current node class
						double sumAccumulator = eln(0);
						for (int Co_cls = 0; Co_cls < cNum; Co_cls++) {
							int max_bitCount = 1 << parentOfCoExcept_currentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label product(fo(p)) except current node
								double productAccumulator = data.allNodes[Co]->fo[Co_cls];
								int parentClsProd = 1;
								for (int p = 0; p < parentOfCoExcept_currentNode.size(); p++) {
									int pid = parentOfCoExcept_currentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									if (pid == Po) {
										double go_Po_child = data.allNodes[Po]->go_ChildList[CIndex * cNum + parentClsValue];
										productAccumulator = elnproduct(productAccumulator, go_Po_child);
									}
									else {
										productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
									}
								}
								//p(Yco|Ypco)
								parentClsProd *= cls;
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[Co_cls][parentClsProd]);
								cout << "Line 848" << endl;
								sumAccumulator = elnsum(sumAccumulator, productAccumulator);
							}
						}
						data.allNodes[cur_node_id]->gi[cls] = sumAccumulator;
					}
				}

				//computing go(n)_parent
				int Co = data.allNodes[cur_node_id]->go_lastVisitedNode;
				for (int cls = 0; cls < cNum; cls++) {
					double productAccumulator = eln(1);
					for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
						int cid = data.allNodes[cur_node_id]->childrenID[c];
						if (cid == Co) {
							continue;
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls]);
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->gi[cls]);
					productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
					data.allNodes[cur_node_id]->go_parent[cls] = productAccumulator;
				}

				//computing go(n)_child
				//for every child c of n . c != Co
				for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
					int cid = data.allNodes[cur_node_id]->childrenID[c];
					if (cid == Co) {
						continue;
					}
					for (int cls = 0; cls < cNum; cls++) {
						double productAccumulator = eln(1);
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi[cls]);
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->gi[cls]);
						for (int d = 0; d < data.allNodes[cur_node_id]->childrenID.size(); d++) {
							if (d == c || data.allNodes[cur_node_id]->childrenID[d] == Co) continue;
							productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[d * cNum + cls]);
						}
						productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
						data.allNodes[cur_node_id]->go_ChildList[c * cNum + cls] = productAccumulator;
					}
				}
			}
			gVisited[cur_node_id] = 1;
			//verification
			for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				if (data.allNodes[cur_node_id]->childrenID[c] == data.allNodes[cur_node_id]->go_lastVisitedNode) {
					for (int cls = 0; cls < cNum; cls++) {
						if (data.allNodes[cur_node_id]->go_ChildList[c * cNum + cls] != 0) {
							cout << " go_childlist Message Computation Error (this should not be computed)  Node " << cur_node_id << endl;
						}
					}
					continue;
				}
				for (int cls = 0; cls < cNum; cls++) {
					if (data.allNodes[cur_node_id]->go_ChildList[c * cNum + cls] < MESSAGELOW || data.allNodes[cur_node_id]->go_ChildList[c * cNum + cls]>0) {
						cout << " go_childlist Message Computation Error in Node " << cur_node_id << endl;
					}
				}
			}
			if (data.allNodes[cur_node_id]->go_fromChild || data.allNodes[cur_node_id]->go_fromParentofChild) {
				for (int cls = 0; cls < cNum; cls++) {
					if (data.allNodes[cur_node_id]->go_parent[cls] < MESSAGELOW || data.allNodes[cur_node_id]->go_parent[cls]>0) {
						cout << " go_parent Message Computation Error in Node " << cur_node_id << endl;
					}
				}
			}
			//verify gi
			for (int cls = 0; cls < cNum; cls++) {
				if (data.allNodes[cur_node_id]->gi[cls] < MESSAGELOW || data.allNodes[cur_node_id]->gi[cls]>0) {
					cout << " gi Message Computation Error in Node " << cur_node_id << endl;
				}
			}
		}
	}
}


//the code assumes 2 by 2 transition matrix P(zn|zpn)
void cFlood::UpdateMarginalProb() {
	// Calculate Marginal distribution

	int curIdx;
	Node* curNode;
	double normFactor;

	for (int i = 0; i < parameter.allPixelSize; i++) {
		curIdx = i;
		curNode = data.allNodes[curIdx];
		//initialize the result variable for cumulation
		for (int zn = 0; zn < cNum; zn++) {  //must initialize
			for (int zpn = 0; zpn < cNum; zpn++) {
				infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn] = LOGZERO;
			}
		}
		// p(z, zp|X, theta) = multiply all outgoing messages towards the factor node (zn|zpn) * p(z|zp)
		// don't forget marginalization over Zpn
		if (curNode->parentsID.size() > 0) {

			if (data.allNodes[curIdx]->go_fromParent) {
				for (int cls = 0; cls < cNum; cls++) {
					int max_bitCount = 1 << data.allNodes[curIdx]->parentsID.size();
					for (int bitCount = 0; bitCount < max_bitCount; bitCount++) {
						double curMessage = data.allNodes[curIdx]->fo[cls];
						int parentClsProd = 1; //p(c), product of parent classes for child c

						for (int p = 0; p < data.allNodes[curIdx]->parentsID.size(); p++) {
							int pid = data.allNodes[curIdx]->parentsID[p];
							int parentClsValue = (bitCount >> p) & 1;
							parentClsProd *= parentClsValue;
							//maintain curMessage with go/fo on parent p
							if (pid == data.allNodes[curIdx]->go_lastVisitedNode) { //Po
								for (int c = 0; c < data.allNodes[pid]->childrenID.size(); c++) {
									if (data.allNodes[pid]->childrenID[c] == curIdx) {
										curMessage = elnproduct(curMessage, data.allNodes[pid]->go_ChildList[c * cNum + parentClsValue]);
										break;
									}
								}
							}
							else {
								curMessage = elnproduct(curMessage, data.allNodes[pid]->fo[parentClsValue]);
							}
						}
						curMessage = elnproduct(curMessage, parameter.elnPz_zpn[cls][parentClsProd]);
						infer.marginal_ZnZpn[curIdx * cNum * cNum + cls * cNum + parentClsProd] = elnsum(infer.marginal_ZnZpn[curIdx * cNum * cNum + cls * cNum + parentClsProd], curMessage);
					}
				}
			}
			else {  //when go is from child or child of parent
				for (int cls = 0; cls < cNum; cls++) {
					//double curMessage = eln(1); //eln(1);
					int max_bitCount = 1 << data.allNodes[curIdx]->parentsID.size();
					for (int bitCount = 0; bitCount < max_bitCount; bitCount++) {
						double curMessage = data.allNodes[curIdx]->go_parent[cls];
						int parentClsProd = 1; //p(c), product of parent classes for child c

						for (int p = 0; p < data.allNodes[curIdx]->parentsID.size(); p++) {
							int pid = data.allNodes[curIdx]->parentsID[p];
							int parentClsValue = (bitCount >> p) & 1;
							parentClsProd *= parentClsValue;
							//maintain curInMessage with go/fo on parent p
							curMessage = elnproduct(curMessage, data.allNodes[pid]->fo[parentClsValue]);
						}
						curMessage = elnproduct(curMessage, parameter.elnPz_zpn[cls][parentClsProd]);
						infer.marginal_ZnZpn[curIdx * cNum * cNum + cls * cNum + parentClsProd] = elnsum(infer.marginal_ZnZpn[curIdx * cNum * cNum + cls * cNum + parentClsProd], curMessage);
					}
				}

			}
			normFactor = LOGZERO;
			for (int zn = 0; zn < cNum; zn++) {
				for (int zpn = 0; zpn < cNum; zpn++) {
					normFactor = elnsum(normFactor, infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn]);
				}
			}

			//marginal_ZnZpn select the first and last term for each z
			for (int zn = 0; zn < cNum; zn++) {
				for (int zpn = 0; zpn < cNum; zpn++) {
					infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn] = elnproduct(infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn], -1 * normFactor);
				}
			}

			//verifying Marginal Probabiltiy sum should be equal to 1
			//value should be in range (-inf,0]
			double sumTest = 0;
			for (int zn = 0; zn < cNum; zn++) {
				for (int zpn = 0; zpn < cNum; zpn++) {
					sumTest += eexp(infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn]);
					if (infer.marginal_ZnZpn[curIdx * cNum * cNum + zn * cNum + zpn] > 0) {
						cout << "Error in Marginal Probability Computation Node " << curIdx << " Zn= " << zn << " Zpn= " << zpn << endl;
					}
				}
			}
			if (abs(sumTest - 1) > 0.0001) {
				cout << "Error in Marginal Probability Computation Node " << curIdx << " Sum is not equal to 1" << endl;
			}

		}

		//else {
		// P(z|X, theta) = gi * fi * vi, Marginal Zn
		normFactor = LOGZERO;
		if (data.allNodes[curIdx]->go_fromParent) {
			for (int z = 0; z < cNum; z++) {
				//compute infer.marginal_Zn[curIdx * cNum + z] based on all incoming messages to the node, including gi's, fi's, and P(xn|Zn)
				//infer.marginal_Zn[curIdx * cNum + z] = elnproduct(elnproduct(infer.lnfi[curIdx * cNum + z], infer.lngi[curIdx * cNum + z]), infer.lnvi[curIdx * cNum + z]);
				double curMessage = data.allNodes[curIdx]->gi[z];
				//incoming from the child side
				for (int c = 0; c < data.allNodes[curIdx]->childrenID.size(); c++) {
					curMessage = elnproduct(curMessage, data.allNodes[curIdx]->fi_ChildList[c * cNum + z]);
				}
				curMessage = elnproduct(curMessage, parameter.elnPxn_zn[curIdx * cNum + z]);
				infer.marginal_Zn[curIdx * cNum + z] = curMessage;
				normFactor = elnsum(normFactor, infer.marginal_Zn[curIdx * cNum + z]);
			}
		}
		else { //when go is from child or child of parent
			for (int z = 0; z < cNum; z++) {
				//compute infer.marginal_Zn[curIdx * cNum + z] based on all incoming messages to the node, including gi's, fi's, and P(xn|Zn)
				//infer.marginal_Zn[curIdx * cNum + z] = elnproduct(elnproduct(infer.lnfi[curIdx * cNum + z], infer.lngi[curIdx * cNum + z]), infer.lnvi[curIdx * cNum + z]);
				double curMessage = data.allNodes[curIdx]->fi[z];
				curMessage = elnproduct(curMessage, data.allNodes[curIdx]->gi[z]);
				//incoming from the child side
				for (int c = 0; c < data.allNodes[curIdx]->childrenID.size(); c++) {
					if (data.allNodes[curIdx]->childrenID[c] == data.allNodes[curIdx]->go_lastVisitedNode) {
						continue;
					}
					curMessage = elnproduct(curMessage, data.allNodes[curIdx]->fi_ChildList[c * cNum + z]);
				}
				curMessage = elnproduct(curMessage, parameter.elnPxn_zn[curIdx * cNum + z]);
				infer.marginal_Zn[curIdx * cNum + z] = curMessage;
				normFactor = elnsum(normFactor, infer.marginal_Zn[curIdx * cNum + z]);
			}
		}
		//}
		for (int z = 0; z < cNum; z++) {
			infer.marginal_Zn[curIdx * cNum + z] = elnproduct(infer.marginal_Zn[curIdx * cNum + z], -1 * normFactor);
		}

		double sumTest = 0;
		for (int c = 0; c < cNum; c++) {
			sumTest += eexp(infer.marginal_Zn[curIdx * cNum + c]);
			if (infer.marginal_Zn[curIdx * cNum + c] > 0) {
				cout << "wrong message: marginal_Zn" << endl;
			}
		}
		if (abs(sumTest - 1) > 0.0001) {
			cout << "Error in Marginal Probability Computation Node " << curIdx << " Sum is not equal to 1" << endl;
		}
	}
}


void cFlood::UpdateParameters() {

	//// Calculate new parameter
	double topEpsilon = LOGZERO, bottomEpsilon = LOGZERO, topPi = LOGZERO, bottomPi = LOGZERO;
	double bottomMu[cNum] = { LOGZERO };
	double tempMu[cNum][Dim] = { LOGZERO };
	double xMinusMu[cNum][Dim];
	double SigmaTemp[cNum][Dim][Dim] = { 0 };

	for (int i = 0; i < parameter.allPixelSize; i++) {
		int curIdx = i;

		// Epsilon, zi has parents
		if (data.allNodes[curIdx]->parentsID.size() > 0) {
			for (int z = 0; z < cNum; z++) {
				for (int zp = 0; zp < cNum; zp++) {
					topEpsilon = elnsum(topEpsilon, elnproduct(eln(zp * (1 - z)), infer.marginal_ZnZpn[i * cNum * cNum + z * cNum + zp]));
					bottomEpsilon = elnsum(bottomEpsilon, elnproduct(eln(zp), infer.marginal_ZnZpn[i * cNum * cNum + z * cNum + zp]));
				}
			}
		}
		// Pi, zi is leaf node
		else {
			for (int z = 0; z < cNum; z++) {
				topPi = elnsum(topPi, elnproduct(eln(1 - z), infer.marginal_Zn[i * cNum + z]));
				bottomPi = elnsum(bottomPi, infer.marginal_Zn[i * cNum + z]);
			}
		}

		// Mu0, Mu1, go through all nodes
		for (size_t j = 0; j < Dim; j++) {
			for (int c = 0; c < cNum; c++) {
				tempMu[c][j] = elnsum(tempMu[c][j], elnproduct(eln(data.features[i * Dim + j]), infer.marginal_Zn[i * cNum + c]));
			}
		}
		for (int c = 0; c < cNum; c++) {
			bottomMu[c] = elnsum(bottomMu[c], infer.marginal_Zn[i * cNum + c]);
		}
	}

	parameter.Epsilon = elnproduct(topEpsilon, -1 * bottomEpsilon);
	parameter.Pi = elnproduct(topPi, -1 * bottomPi);


	// reserve eln(Mu) form
	for (size_t j = 0; j < Dim; j++) {
		for (int c = 0; c < cNum; c++) {
			parameter.elnMu[c][j] = elnproduct(tempMu[c][j], -1 * bottomMu[c]);
		}
	}

	// convert Mu to normal
	for (size_t j = 0; j < Dim; j++) {
		for (int c = 0; c < cNum; c++) {
			parameter.Mu[c][j] = eexp(parameter.elnMu[c][j]);
		}
	}


	// Update Sigma
	for (size_t i = 0; i < parameter.allPixelSize; i++) {

		for (int c = 0; c < cNum; c++) {
			for (size_t j = 0; j < Dim; j++) {
				xMinusMu[c][j] = data.features[i * Dim + j] - parameter.Mu[c][j];
			}
		}

		for (int c = 0; c < cNum; c++) {
			for (size_t m = 0; m < Dim; m++) { // row
				for (size_t n = 0; n < Dim; n++) { // column
					SigmaTemp[c][m][n] += xMinusMu[c][m] * xMinusMu[c][n] * eexp(infer.marginal_Zn[i * cNum + c]);
				}
			}
		}
	}

	for (int c = 0; c < cNum; c++) {
		for (size_t i = 0; i < Dim; i++) {
			for (size_t j = 0; j < Dim; j++) {
				parameter.Sigma[c][i][j] = SigmaTemp[c][i][j] / eexp(bottomMu[c]); // bottom is the same as Mu
			}
		}
	}

}

vector<int> cFlood::getBFSOrder(int root, vector<int>& bfsVisited, int bank) {
	//vector<int> bfsVisited;
	vector<int> bfs;
	queue<int> que;
	que.push(root);

	while (!que.empty()) {
		int currentNode = que.front();
		bfs.push_back(currentNode);
		bfsVisited[currentNode] = 1;
		que.pop();
		for (int i = 0; i < data.allNodes[currentNode]->childrenID.size(); i++) {
			int child = data.allNodes[currentNode]->childrenID[i];
			if (!bfsVisited[child]) {
				//if (data.allNodes[child]->bank == 0 || data.allNodes[child]->bank == bank) {
				que.push(child);
				//}

			}
		}
		for (int i = 0; i < data.allNodes[currentNode]->parentsID.size(); i++) {
			int parent = data.allNodes[currentNode]->parentsID[i];
			if (!bfsVisited[parent]) {
				//if (data.allNodes[parent]->bank == 0 || data.allNodes[parent]->bank == bank) {
				que.push(parent);
				//}
				//que.push(parent);
			}
		}
	}
	return bfs;
}



void cFlood::input(int argc, char* argv[]) {
	cout << "argc: " << argc << endl;

	GDALAllRegister();
	clock_t start_s = clock();

	cout << "start here " << endl;

	if (argc > 1) {

		ifstream config(argv[1]);
		string line;
		getline(config, line);
		CTInputLocation = line;  //Input file location
		getline(config, line);
		CTSourceDirection = line;           //Elevation data file name
		getline(config, line);
		CTProbability = line;
		//getline(config, line);
		//CTBank = line;
		//getline(config, line);
		//CTCost = line;
		//getline(config, line);
		//CTPits = line;
		//getline(config, line);
		//CTTree = line;
		//getline(config, line);
		//CTRoughness = line;
		getline(config, line);
		CTFel = line;
		getline(config, line);
		CTPara = line;          //parameter data file name
		getline(config, line);
		CTStream = line; // raster file for North Carolina reach
		getline(config, line);
		CTOutputFolderByDate = line; //oputput folder by date
		getline(config, line);
		CTOutputLocation = line; //oputput location to store the output of HMCT
		getline(config, line);
		CTPrediction = line;    //file name for output prediction data
// 		getline(config, line);
// 		CTLeftBank = line;
		getline(config, line);
		CTPredictionTxt = line;
		getline(config, line);
		CTPredictionRegularized = line;
		getline(config, line);
		CTPredictionRegularizedTxt = line;
		getline(config, line);
		CTTestCase = line;
		getline(config, line);
		CTAlgo = line;


	}
	else {
		std::cout << "Missing Configuration File!";
	}

	struct stat info;

	cout << "CTInputLocation: " << CTInputLocation << endl;


	int status = dirExists(CTInputLocation.c_str());
	cout << "status: " << status << endl;
	if (status <= 0) {
		cout << "Error: input directory does not exist.." << endl;
		exit(0);
	}

	status = dirExists(CTOutputFolderByDate.c_str());
	if (status <= 0) {
		cout << "Error: output folder by date does not exist..creating one!" << endl;
		//exit(0);
		status = mkdir(CTOutputFolderByDate.c_str(), 0777); // Added by Saugat: create output dir if not exists
		if (status != 0) {
			cout << "Error: could not create Output Folder by date.." << endl;
			exit(0);
		}
	}


	status = dirExists(CTOutputLocation.c_str());
	cout << " Output Location: " << CTOutputLocation << endl;
	if (status <= 0) {
		cout << "Error: output directory does not exist..creating one!" << endl;
		//exit(0);
		status = mkdir(CTOutputLocation.c_str(), 0777); // Added by Saugat: create output dir if not exists
		if (status != 0) {
			cout << "Error: could not create Output Directory.." << endl;
			exit(0);
		}
	}
	status = dirExists((CTOutputLocation + "ProfileTables").c_str());
	if (status <= 0) {
		status = mkdir((CTOutputLocation + "ProfileTables").c_str(), 0777);
		if (status != 0) {
			cout << "Error: could not create ProfileTables folder.." << endl;
			exit(0);
		}
		//exit(0);
	}

	status = dirExists((CTOutputLocation + "BoundaryTables").c_str());
	if (status <= 0) {
		status = mkdir((CTOutputLocation + "BoundaryTables").c_str(), 0777);
		if (status != 0) {
			cout << "Error: could not create BoundaryTables folder.." << endl;
			exit(0);
		}
	}




	//reading text file
	ifstream parameterFile(CTInputLocation + CTPara);
	if (!parameterFile) {
		std::cout << "Failed to open parameter!" << endl;
		exit(0);
	}
	parameterFile >> parameter.reachId;
	parameterFile >> parameter.Epsilon;
	parameterFile >> parameter.Pi;
	parameterFile >> parameter.cutoff_percentage;
	parameterFile >> parameter.minCost;
	parameterFile >> parameter.useCutoff;
	parameterFile >> parameter.useHMT;
	parameterFile >> parameter.lambda;
	parameterFile >> parameter.num_range;
	parameter.Pi_orig = parameter.Pi;
	parameter.lambda_str.precision(2);
	parameter.lambda_str << std::fixed << parameter.lambda;

	//parameter.reachId = CTPara.substr(9, 2);
	parameterFile.close();



	GeotiffRead sourceDirTiff((CTInputLocation + CTSourceDirection).c_str());
	cout << "after src dir" << endl;
	GeotiffRead floodProbTiff((CTInputLocation + CTProbability).c_str());
	GeotiffRead reachTiff((CTInputLocation + CTStream).c_str());

	GeotiffRead FelTiff((CTInputLocation + CTFel).c_str());


	if (parameter.useHMT) {
		cout << "using HMT tree" << endl;
	}

	// The pointer to the raster band data of the source direction tiff
	float** sourceDirData = (float**)sourceDirTiff.GetRasterBand(1);
	float** floodProbData = floodProbTiff.GetRasterBand(1);

	float** felData = FelTiff.GetRasterBand(1);
	float** reachData = reachTiff.GetRasterBand(1);


	// Get the array dimensions
	int* dims = sourceDirTiff.GetDimensions();
	// Get the nodata value
	//float noDataValue = (double)sourceDirTiff.GetNoDataValue();

	parameter.ROW = dims[0];
	parameter.COLUMN = dims[1];
	parameter.allPixelSize = parameter.ROW * parameter.COLUMN;
	parameter.orgPixelSize = parameter.allPixelSize;

	//Note end

	if (parameter.Epsilon > 1 || parameter.Pi > 1) {
		std::cout << "wrong parameter" << endl;
	}

	std::cout << "Input parameters:" << endl << "Reach Id:" << parameter.reachId << "Epsilon: " << parameter.Epsilon << " Pi: " << parameter.Pi << endl;
	cout << "rho: " << parameter.rho << endl;


	// TODO: added by Saugat
	data.features.resize(parameter.allPixelSize * Dim);// RGB + ..., rowwise, long array
	parameter.NAValue = -1; // TODO: check

	data.NA.resize(parameter.allPixelSize);
	std::fill(data.NA.begin(), data.NA.end(), false);
	for (size_t i = 0; i < parameter.allPixelSize; i++) {
		for (int j = 0; j < Dim; j++) {
			if (data.features[i * Dim + j] == parameter.NAValue) {
				data.NA[i] = true;
				break;
			}
		}
	}

	// Fill in the data values in the tree
	int nodeIdCounter = 0;
	// 	std::map<int, int> nodeIndexMap;
	int NCOLS = dims[1];
	int NROWS = dims[0];

	int skipped_nodes = 0;

	cout << "NROWS: " << NROWS << " NCOLS: " << NCOLS << endl;

	//data.allNodes.resize(NROWS * NCOLS);
	int src_dir_m_1 = 0;
	int tmp_idx = 0;
	int tmp_idx_2 = 0;
	for (int row = 0; row < NROWS; row++)
	{
		for (int col = 0; col < NCOLS; col++)
		{
			////cout << tmp_idx;
			//if (tmp_idx % 1000 == 0) {
			//	cout << tmp_idx << "  ";
			//}

			int sourceDir = sourceDirData[row][col];
			float floodProb = floodProbData[row][col];
			//int bank = bankData[row][col];
			// int cost = costData[row][col]; // TODO: check

			//double cost = costData[row][col];
			//int pits = (pitsData[row][col] == 1) ? 1 : 0;
			//int tree = (treeData[row][col] > 0) ? 1 : 0;
			//int roughness = (roughnessData[row][col] >= 0.5) ? 1 : 0;
			float fel = felData[row][col];

			int rch = reachData[row][col];

			if (rch > 0) { // if -1 is not set, put >0 else >= 0
				src_dir_m_1++;
				sourceDir = 0;
			}

			int currentId = row * NCOLS + col;
			int parentId;
			bool river_node = false;

			// get left bank data
// 			int isLeftBankNode = leftBankData[row][col];
// 			if (isLeftBankNode == -11) {
// 				data.leftBankNodes.insert(make_pair(currentId, true));
// 			}
// 		/*	if (row == 0 || row == NROWS - 1 || col == 0 || col == NCOLS - 1)
// 			{
// 				sourceDir = -1;
// 			}*/
			switch (sourceDir) {
			case 1:
				parentId = currentId + 1;
				break;
			case 8:
				parentId = currentId + parameter.COLUMN + 1;
				break;
			case 7:
				parentId = currentId + parameter.COLUMN;
				break;
			case 6:
				parentId = currentId + parameter.COLUMN - 1;
				break;
			case 5:
				parentId = currentId - 1;
				break;
			case 4:
				parentId = currentId - parameter.COLUMN - 1;
				break;
			case 3:
				parentId = currentId - parameter.COLUMN;
				break;
			case 2:
				parentId = currentId - parameter.COLUMN + 1;
				break;
			case 0:
				data.reaches.push_back(currentId);
				data.reaches_map.insert(make_pair(currentId, true));
				river_node = true;
				break;
				//continue;
			default:
				skipped_nodes++; // only no data cases // only boundary nodes have no data // middle nodes have both parent and child
				// continue; // issue with flow direction layer
				break;
			}
			int newCurrentId, newParentId;
			if (nodeIndexMap.find(currentId) == nodeIndexMap.end()) {
				nodeIndexMap.insert(make_pair(currentId, nodeIdCounter));
				newCurrentId = nodeIdCounter;
				data.allNodes.push_back(new Node(0, nodeIdCounter));
				data.allNodes[newCurrentId]->originalId = currentId;
				nodeIdCounter++;
				tmp_idx_2++;

			}
			else {
				newCurrentId = nodeIndexMap[currentId];
				tmp_idx_2++;
			}



			// TODO: check: River ids with no parents
			if (river_node == false) {
				if (nodeIndexMap.find(parentId) == nodeIndexMap.end()) {
					nodeIndexMap.insert(make_pair(parentId, nodeIdCounter));
					//nodeIndexMap[parentId] = nodeIdCounter;
					newParentId = nodeIdCounter;
					data.allNodes.push_back(new Node(0, nodeIdCounter));
					data.allNodes[newParentId]->originalId = parentId;
					nodeIdCounter++;
				}
				else {
					newParentId = nodeIndexMap[parentId];
				}

				data.allNodes[newCurrentId]->parentsID.push_back(newParentId);   //source direction layer
				data.allNodes[newParentId]->childrenID.push_back(newCurrentId);  //source direction layer
			}

			data.allNodes[newCurrentId]->isObserved = (floodProb == 0) ? 0 : 1;           // no pits no na

			data.allNodes[newCurrentId]->isNa = (floodProb == 0) ? 1 : 0;                 //na?  0-255 deltaresult.tif  0 is na

			data.allNodes[newCurrentId]->p = floodProb / 255.0;                   // delta
			data.allNodes[newCurrentId]->fel = fel;                 //field elevation layer

			tmp_idx++;

		}
	}
	cout << "Start Debugging" << endl;

	cout << nodeIndexMap.size() << endl;
	cout << "counter: " << data.allNodes.size() << endl;
	cout << "row cols total: " << NROWS * NCOLS << endl;

	cout << "skipped: " << skipped_nodes << endl;
	cout << "counter: " << data.allNodes.size() << endl;
	cout << "NROWS: " << NROWS << " NCOLS: " << NCOLS << endl;
	cout << "Src dir -1: " << src_dir_m_1 << endl;
	cout << "reaches size: " << data.reaches.size() << endl;
	cout << "tmp idx: " << tmp_idx << endl;
	cout << "tmp idx_2: " << tmp_idx_2 << endl;
	cout << "node id counter: " << nodeIdCounter << endl;

	cout << "total size: " << data.allNodes.size() << endl;



	// get the reach nodes
	for (int i = 0; i < data.reaches.size(); i++) {
		int currentreachId = data.reaches[i];
		int newReachId = nodeIndexMap[currentreachId];

		if (data.allNodes[newReachId]->parentsID.size() == 0 && data.allNodes[newReachId]->childrenID.size() != 0) {
			data.reach_ids.push_back(newReachId);
			data.reach_ids_orig.push_back(currentreachId);

			data.reach_ids_map.insert(make_pair(newReachId, true));
			data.reach_ids_orig_map.insert(make_pair(currentreachId, true));

			////
			int row = (int)(currentreachId / parameter.COLUMN);
			int col = currentreachId % parameter.COLUMN;

			float fel = felData[row][col];
			data.reach_fel.push_back(fel);
			int newCurrentId = -1;
			if (nodeIndexMap.find(currentreachId) == nodeIndexMap.end()) {
				cout << "nodeIdCounter: " << nodeIdCounter << endl;
				nodeIndexMap.insert(make_pair(currentreachId, nodeIdCounter));
				//nodeIndexMap[currentId] = nodeIdCounter;
				newCurrentId = nodeIdCounter;
				data.allNodes.push_back(new Node(0, nodeIdCounter));
				data.allNodes[newCurrentId]->originalId = currentreachId;
				data.AdjustedReachNodes.push_back(newCurrentId);
				// data.allNodes[newCurrentId]->bank = 0;
				data.allNodes[newCurrentId]->cost = 0.0;
				data.allNodes[newCurrentId]->p = 1.0; // check
				data.allNodes[newCurrentId]->isObserved = 0; // to indicate reach nodes
				data.allNodes[newCurrentId]->fel = fel;
				nodeIdCounter++;

			}
			else {
				newCurrentId = nodeIndexMap[currentreachId];
				data.AdjustedReachNodes.push_back(newCurrentId);
				// data.allNodes[newCurrentId]->bank = 0;
				data.allNodes[newCurrentId]->cost = 0.0;
				data.allNodes[newCurrentId]->p = 1.0; // check
				data.allNodes[newCurrentId]->isObserved = 0; // to indicate reach nodes
				data.allNodes[newCurrentId]->fel = fel;
			}
			////

		}
		else if (data.allNodes[newReachId]->parentsID.size() == 0 && data.allNodes[newReachId]->childrenID.size() == 0) {
			data.river_ids.push_back(newReachId);
			data.river_ids_map.insert(make_pair(currentreachId, true));
		}
		else {
			data.river_ids_orig.push_back(newReachId);
		}
	}

	cout << "AdjustedReachNodes: " << data.AdjustedReachNodes.size() << endl;
	cout << "reach ids: " << data.reach_ids.size() << endl;



	cout << "river size: " << data.river_ids.size() << endl;
	cout << "reaches size: " << data.reach_ids.size() << endl;
	cout << "river ids orig: " << data.river_ids_orig.size() << endl;

	// TODO: 1. check tree
	// 2. Find reaches based on parent-child relation and save reach ids
	// 3. separate regions

	//adjust reach ids with new index structure
	extra.reachLength = 0;
	double rowSum = 0;
	double colSum = 0;



	// CHECK: get left and right bank source nodes in order of stream
	std::map<int, bool> bfs_visited;
	int tmp_node_id;
	queue<pair<int, int>> bfs_que;
	std::map<int, bool> on_queue;
	vector<int> node_order;
	vector<int> left_node_order;
	vector<int> right_node_order;


	for (int row = 0; row < parameter.ROW; row++) {
		for (int col = 0; col < parameter.COLUMN; col++) {
			int node_id_curr = row * parameter.COLUMN + col;
			// skip nodes not in the river
			if (!data.reach_ids_orig_map[node_id_curr] && !data.river_ids_map[node_id_curr]) {
				continue;
			}


			tmp_node_id = row * parameter.COLUMN + col;

			if (bfs_visited.find(tmp_node_id) != bfs_visited.end()) {
				continue;
			}

			bfs_que.push(make_pair(row, col));
			bfs_visited[tmp_node_id] = true;
			on_queue[tmp_node_id] = true;

			while (!bfs_que.empty()) {
				pair<int, int> curr_node = bfs_que.front();
				int i = curr_node.first;
				int j = curr_node.second;

				int idx = i * parameter.COLUMN + j;

				bfs_visited[idx] = true;
				node_order.push_back(idx);
				bfs_que.pop();

				int l, r;
				for (l = -1; l < 2; l++) {
					for (r = -1; r < 2; r++) {
						if (l == 0 && r == 0) {
							continue;
						}

						int i_nei = i + l;
						int j_nei = j + r; // get the neighboring x and y

						// check for boundary cases
						if (i_nei < 0 || j_nei < 0 || j_nei >= parameter.COLUMN || i_nei >= parameter.ROW) {
							continue;
						}

						// check if already visited or not
						int neigh_node_id = i_nei * parameter.COLUMN + j_nei;

						if (bfs_visited.find(neigh_node_id) != bfs_visited.end()) {
							continue;
						}

						if (on_queue.find(neigh_node_id) != on_queue.end()) {
							continue;
						}

						if (data.reach_ids_orig_map[neigh_node_id] || data.river_ids_map[neigh_node_id]) {
							on_queue[neigh_node_id] = true;
							bfs_que.push(make_pair(i_nei, j_nei));
						}

					}

				}
			}

		}
	}





	// get left and right


	for (int i = 0; i < node_order.size(); i++) {
		int tmp_node_ = node_order[i];
		if (data.river_ids_map[tmp_node_]) {
			continue;
		}

		data.rightNodesInOrder.push_back(nodeIndexMap[tmp_node_]);


	}

	cout << "Node order size: " << node_order.size() << endl;

	cout << "Right nodes size: " << data.rightNodesInOrder.size() << endl;


	data.inferedmaxCostRight.resize(data.rightNodesInOrder.size(), -1);


	data.regularizedMaxCostRight.resize(data.rightNodesInOrder.size(), -1);



	data.rightbfsRootNodes.resize(data.rightNodesInOrder.size(), -1); //leaf nodes with no children


	// Right
	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		int nodeId = data.rightNodesInOrder[i];
		//finding root in each tree
		int rightnode = nodeId;



		data.rightbfsRootNodes[i] = rightnode;
	}
	//get bfs order for each tree



	vector<int> bfsVisited;
	bfsVisited.resize(data.allNodes.size(), 0);



	// Right
	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		if (data.rightbfsRootNodes[i] == -1) {
			data.rightbfsOrder.push_back({});
		}
		else {
			data.rightbfsOrder.push_back(getBFSOrder(data.rightbfsRootNodes[i], bfsVisited, 2));
		}
	}

	// Check size 
	cout << "All Nodes in Order: " << endl;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		int rgn_id = data.rightNodesInOrder[rightOrder];
		if (data.rightbfsOrder[rightOrder].size() < PIXELLIMT){
			continue;
		}

		cout << rgn_id << endl;
	}


	// Commented out writing region map

	// write pixels -> region id to raster file
	cout << "Writing region map to tiff" << endl;
	vector<int> map_region_id;
	map_region_id.resize(parameter.orgPixelSize, -1);


	int px_id;
	int rgn_id;


	// Right
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		rgn_id = data.rightNodesInOrder[rightOrder];

		double regionRowSum = 0;
		double regionColSum = 0;
		for (int p = 0; p < data.rightbfsOrder[rightOrder].size(); p++) {
			px_id = data.rightbfsOrder[rightOrder][p];
			px_id = data.allNodes[px_id]->originalId;
			map_region_id[px_id] = rgn_id;


		}


	}



	// TODO: remove
	int minus_2_count = 0;
	int one_count = 0;

	float** region_map = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		region_map[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{

			if (data.reach_ids_orig_map[index]) {
				// use -2 for boundary nodes in middle of river
				minus_2_count++;
				region_map[row][col] = -2;
			}
			else if (data.river_ids_map[index]) {
				// use 1 for nodes in river
				one_count++;
				region_map[row][col] = 1;
			}
			else {
				region_map[row][col] = map_region_id[index];

			}
			index++;
		}
	}

// 	cout << "minus_2_count: " << minus_2_count << endl;
// 	cout << "one_count: " << one_count << endl;



// 	cout << "Reach ids: " << data.reach_ids_orig.size() << endl;

	GDALDataset* srcDataset_ = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform_[6];
	srcDataset_->GetGeoTransform(geotransform_);
	const OGRSpatialReference* poSRS_ = srcDataset_->GetSpatialRef();

	GeotiffWrite mapTiff((CTOutputLocation + parameter.reachId + "_RegionMap" + ".tif").c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform_, poSRS_);
	mapTiff.write(region_map);


	cout << "Writing region map to tiff completed!!!" << endl;
	// region map comment end

// 	return;

	// data.hasObservedPixelsLeft.resize(data.leftbfsOrder.size(), 0);
	data.hasObservedPixelsRight.resize(data.rightbfsOrder.size(), 0);



	// Right
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->isObserved == 1) {
				data.hasObservedPixelsRight[i] = 1;
				break;
			}
		}
	}

	vector<float> cost_map;
	cost_map.resize(parameter.orgPixelSize, -1);



	// Right
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		int reach_id_current = data.rightNodesInOrder[i];
		float reach_fel = data.allNodes[reach_id_current]->fel;
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			float current_fel = data.allNodes[pixelId]->fel;
			data.allNodes[pixelId]->cost = current_fel - reach_fel;
			cost_map[data.allNodes[pixelId]->originalId] = current_fel - reach_fel;
		}
	}

	float** cost_data = new float* [parameter.ROW];
	int indexx = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		cost_data[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
			cost_data[row][col] = cost_map[indexx];
			indexx++;
		}
	}

	cout << "Writing cost map to tiff!!!" << endl;
	// save cost as a new tif file
	GDALDataset* srcDataset_3 = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform_3[6];
	srcDataset_3->GetGeoTransform(geotransform_3);
	const OGRSpatialReference* poSRS_3 = srcDataset_3->GetSpatialRef();

	GeotiffWrite mapTiff3((CTOutputLocation + parameter.reachId + "_CostMap" + ".tif").c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform_3, poSRS_3);
	mapTiff3.write(cost_data);
	cout << "Writing cost map to tiff completed!!!" << endl;




	////// TODO: call function to export
	////export_FIST_structure("Child_list_Reach19_FIST.txt", "Parent_list_Reach19_FIST.txt", "Small_Region_Pixels_FIST.txt");


	// data.highestCostLeft.resize(data.leftNodesInOrder.size(), MAXCOST);
	data.highestCostRight.resize(data.rightNodesInOrder.size(), MAXCOST);



	// Right
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->cost > data.highestCostRight[i]) {
				data.highestCostRight[i] = data.allNodes[pixelId]->cost;
			}
		}
	}
	parameter.allPixelSize = data.allNodes.size();
	parameter.elnPzn_xn.resize(parameter.allPixelSize * cNum, eln(0.5));




	//// sort for trees on the right bank
	data.costIndexPairRight.resize(data.rightbfsOrder.size());
	data.sortedCostIndexRight.resize(data.rightbfsOrder.size());
	data.rightIndex2PixelId.resize(data.rightbfsOrder.size());

	for (int i = 0; i < data.rightbfsOrder.size(); i++) {

		data.costIndexPairRight[i].resize(data.rightbfsOrder[i].size());
		data.sortedCostIndexRight[i].resize(data.rightbfsOrder[i].size());

		for (size_t j = 0; j < data.rightbfsOrder[i].size(); j++) {
			int pixelId = data.rightbfsOrder[i][j];
			data.rightIndex2PixelId[i][j] = pixelId;
			data.costIndexPairRight[i][j] = make_pair(data.allNodes[pixelId]->cost, j); // sort by cost
			//data.costIndexPairRight[i][j] = make_pair(data.allNodes[pixelId]->fel, j); // sort by elevation
		}
		sort(std::begin(data.costIndexPairRight[i]), std::end(data.costIndexPairRight[i]));

		for (int k = 0; k < data.rightbfsOrder[i].size(); k++) {
			int index = data.costIndexPairRight[i][k].second;
			data.sortedCostIndexRight[i][index] = k;
		}
	}





	if (parameter.useHMT) {
		for (int i = 0; i < parameter.allPixelSize; i++) {
			data.allNodes[i]->childrenID.clear();
			data.allNodes[i]->parentsID.clear();
		}


		splitTree();


	}

	parameter.elnPzn_xn.resize(parameter.allPixelSize * cNum, eln(0.5));
	for (int i = 0; i < parameter.allPixelSize; i++) {
		if (data.allNodes[i]->isObserved == -1) {
			parameter.elnPzn_xn[i * cNum] = eln(0.5);
			parameter.elnPzn_xn[i * cNum + 1] = eln(0.5);
		}
		else if (data.allNodes[i]->isObserved == 0) {  //stream nodes

			// Old code
			data.allNodes[i]->p = 0.999;
			parameter.elnPzn_xn[i * cNum] = eln(0.001);
			parameter.elnPzn_xn[i * cNum + 1] = eln(0.999);
		}
		else {
			float probability = data.allNodes[i]->p;
			if (probability < 0 || probability>1) {
				std::cout << "i= " << i << " prob erorr :" << probability << endl;
			}
			parameter.elnPzn_xn[i * cNum] = eln(1 - probability);
			parameter.elnPzn_xn[i * cNum + 1] = eln(probability);
		}

		// TODO: test for isNA regions issue; isNA should have low p value
		if (data.allNodes[i]->isNa == 1) {  //NA nodes
			data.allNodes[i]->p = 0.5;
			parameter.elnPzn_xn[i * cNum] = eln(0.5);
			parameter.elnPzn_xn[i * cNum + 1] = eln(0.5);
		}
	}


	parameter.Pi = eln(parameter.Pi);
	parameter.Epsilon = eln(parameter.Epsilon); //check if already eln form?



	cout << "before update transprob" << endl;

	this->UpdateTransProb();
	string reachId = CTPara.substr(0, 2);
	parameter.reachId = reachId;
	string EffectiveBrach = "CL";
	if (EB == 1) {
		EffectiveBrach = "MC";
	}
	else if (EB == 2) {
		EffectiveBrach = "BC";
	}
	parameter.fname = reachId + "_pr_DV3_Top" + to_string((int)(parameter.cutoff_percentage * 100)) + "_" + EffectiveBrach + "_BN" + to_string(BOUNDARY_NODES_OBSERVED);



	std::cout << "Inference Started.." << endl;
	auto start = std::chrono::system_clock::now();
	inference();



	// 	map<int, bool> large_regions;
	map<int, int> id2idx;


	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		id2idx.insert(make_pair(data.rightNodesInOrder[i], i));
		if (data.inferedmaxCostRight[i] > -1) {
			data.large_regions.insert(make_pair(data.rightNodesInOrder[i], true));
		}
	}

	cout << "Large Regions size: " << data.large_regions.size() << endl;
	cout << "Left Nodes in Order: " << endl;



// 	cout << "Right Nodes in Order: " << endl;
	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		// pair<int, int> it;
		 // int it = id2idx[data.rightNodesInOrder[i]];
		double cost_ = data.inferedmaxCostRight[i];
		// data.inferedmaxCostRight_new.push_back(cost_);
		data.inferedmaxCost_id2cost.insert(make_pair(data.rightNodesInOrder[i], cost_));
		// data.inferedmaxCostRight_idx2id.insert(make_pair(i, data.rightNodesInOrder[i]));

		if (data.large_regions[data.rightNodesInOrder[i]]) {
// 			cout << data.rightNodesInOrder[i] << endl;
			data.rightNodesInCorrectOrder.push_back(data.rightNodesInOrder[i]);
		}
	}



	cout << "before selected prediction fist" << endl;

	//////delta_prediction();

	if (parameter.useHMT) {
		selected_prediction();
	}
	else {
		selected_prediction_FIST();
	}

	cout << "after selected prediction fist" << endl;

	interpolate();

	cout << "after interpolate" << endl;

	// we don't need to predict again for HMT
	if (parameter.useHMT) {
		prediction();
	}
	else {
		prediction_FIST();
	}
	output();


	// TODO: check, calculate loglikelihood for regularization
	getLoglikelihood();
	vector<int> regularizedSizeLeft, regularizedSizeRight;
	ofstream outputfile;
	// for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++)
	// {
	// 	if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
	// 		continue;
	// 	}
	// 	regularizedSizeLeft.push_back(data.leftbfsOrder[leftOrder].size());
	// }
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++)
	{
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		regularizedSizeRight.push_back(data.rightbfsOrder[rightOrder].size());
	}

	
	timeSave.open("timeIntervals.txt");

    viterbi(parameter.lambda,parameter.num_range);
		

	outputRegulization(std::to_string(parameter.lambda));
	////cout << "Range Value Now" << ranges[i] << endl;
	data.regularizedMaxCostRight.clear();
	data.regularizedMaxCostRight.resize(data.rightNodesInOrder.size(), -1);
	

}



int cFlood::reachBFS(queue<pair<int, int>>& bfs_que, map<int, bool>& bfs_visited, vector<int>& left_node_order, map<int, bool>& on_queue) {
	int last_node;
	while (!bfs_que.empty()) {
		pair<int, int> curr_node = bfs_que.front();

		int i = curr_node.first;
		int j = curr_node.second;

		int idx = i * parameter.COLUMN + j;
		last_node = idx;
		bfs_visited.insert(make_pair(idx, true));

		if (data.reach_ids_orig_map[idx]) {
			//  cout << "leftttt: " << idx << endl;
			left_node_order.push_back(idx);
		}

		bfs_que.pop();

		for (int l = -1; l < 2; l++) {
			for (int r = -1; r < 2; r++) {
				if (l == 0 && r == 0) {
					continue;
				}

				int i_nei, j_nei = (i + l, j + r); // get the neighboring x and y

				// check for boundary cases
				if (i_nei < 0 || j_nei < 0 || j_nei >= parameter.COLUMN || i_nei >= parameter.ROW) {
					continue;
				}

				// check if already visited or not
				int neigh_node_id = i_nei * parameter.COLUMN + j_nei;

				// skip nodes not on left side
				// if (!data.leftBankNodes[neigh_node_id]) {
				// 	continue;
				// }

				if (bfs_visited[neigh_node_id]) {
					continue;
				}

				if ((data.reach_ids_orig_map[neigh_node_id]) && (!on_queue[neigh_node_id])) {
					bfs_que.push(make_pair(i_nei, j_nei));
					on_queue.insert(make_pair(neigh_node_id, true));
				}
			}
		}
	}

	return last_node;
}




void cFlood::export_FIST_structure(string child_file, string parent_file, string small_file) {
	ofstream cfile;
	cfile.open(child_file);

	ofstream pfile;
	pfile.open(parent_file);

	ofstream sfile;
	sfile.open(small_file);

	cout << "Starting Left\n";

	int Leftcount = 0;
	for (int i = 0; i < data.leftbfsOrder.size(); i++) {

		// skip small regions
		if (data.leftbfsOrder[i].size() < PIXELLIMT) {
			// add pixels of small region to file: Jiaqing needs this
			for (int idx = 0; idx < data.leftbfsOrder[i].size(); idx++) {
				int pxlId = data.leftbfsOrder[i][idx];
				sfile << pxlId;
				sfile << endl;
			}

			continue;
		}


		Leftcount = Leftcount + data.leftbfsOrder[i].size();
		for (int treeIndex = 0; treeIndex < data.leftbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.leftbfsOrder[i][treeIndex];
			cfile << pixelId;
			pfile << pixelId;
			if (data.allNodes[pixelId]->childrenID.size() > 0)
				cfile << ",";
			if (data.allNodes[pixelId]->parentsID.size() > 0)
				pfile << ",";

			for (int k = 0; k < data.allNodes[pixelId]->childrenID.size(); k++) {
				if (k == (data.allNodes[pixelId]->childrenID.size() - 1))
					cfile << data.allNodes[pixelId]->childrenID[k];
				else
					cfile << data.allNodes[pixelId]->childrenID[k] << ",";
			}


			for (int k = 0; k < data.allNodes[pixelId]->parentsID.size(); k++) {
				if (k == (data.allNodes[pixelId]->parentsID.size() - 1))
					pfile << data.allNodes[pixelId]->parentsID[k];
				else
					pfile << data.allNodes[pixelId]->parentsID[k] << ",";
			}

			cfile << endl;
			pfile << endl;
		}
	}

	int Rightcount = 0;
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {

		// skip small regions
		if (data.rightbfsOrder[i].size() < PIXELLIMT) {
			// add pixels of small region to file: Jiaqing needs this
			for (int idx = 0; idx < data.rightbfsOrder[i].size(); idx++) {
				int pxlId = data.rightbfsOrder[i][idx];
				sfile << pxlId;
				sfile << endl;
			}

			continue;
		}


		Rightcount = Rightcount + data.rightbfsOrder[i].size();
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			cfile << pixelId;
			pfile << pixelId;
			if (data.allNodes[pixelId]->childrenID.size() > 0) {
				cfile << ",";
			}
			if (data.allNodes[pixelId]->parentsID.size() > 0) {
				pfile << ",";
			}

			for (int k = 0; k < data.allNodes[pixelId]->childrenID.size(); k++) {
				if (k == (data.allNodes[pixelId]->childrenID.size() - 1))
					cfile << data.allNodes[pixelId]->childrenID[k];
				else
					cfile << data.allNodes[pixelId]->childrenID[k] << ",";
			}

			for (int k = 0; k < data.allNodes[pixelId]->parentsID.size(); k++) {
				if (k == (data.allNodes[pixelId]->parentsID.size() - 1))
					pfile << data.allNodes[pixelId]->parentsID[k];
				else
					pfile << data.allNodes[pixelId]->parentsID[k] << ",";
			}

			cfile << endl;
			pfile << endl;
		}
	}

	cout << Rightcount << " " << Leftcount << endl;
	cout << "Done\n";

	cfile.close();
	pfile.close();

	return;
}

// function to do prediction only on selected regions for easy validation
// commented the water filling part, don't predict again, just use the inferred class
void cFlood::selected_prediction() {
	cout << "Selected prediction started!" << endl;

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		/*if (data.hasObservedPixelsRight[rightOrder] && data.rightbfsOrder[rightOrder].size()>= PIXELLIMT) {
			continue;
		}*/
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];

				// commenting temporarily
				if (data.inferedmaxCostRight[rightOrder] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					continue;
					/*if (data.allNodes[nodid]->isNa == 0)
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder] * 2;*/
				}
				else {

					// Refill
					if (data.allNodes[nodid]->cost <= data.inferedmaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) { // CHECK NC
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder];
					}
					else {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}




				}
			}
		}
	}

	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		//mappredictions[data.rightNodesInOrder[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.rightNodesInOrder[i]]->originalId] = 1;
	}

	// nodes in river(high prob pixels from UNet should also be flooded)
	for (int i = 0; i < data.river_ids.size(); i++) {
		//mappredictions[data.river_ids_orig[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.river_ids[i]]->originalId] = 1;
	}


}

// function to do prediction only on selected regions for easy validation
// commented the water filling part, don't predict again, just use the inferred class
void cFlood::selected_prediction_FIST() {
	cout << "Selected prediction started!" << endl;


	for (int i = 0; i < data.reach_ids.size(); i++) {
		//mappredictions[data.reach_ids[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.reach_ids[i]]->originalId] = 1;
	}

	// nodes in river(high prob pixels from UNet should also be flooded)
	for (int i = 0; i < data.river_ids.size(); i++) {
		//mappredictions[data.river_ids_orig[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.river_ids[i]]->originalId] = 1;
	}

	// Comment: END

	auto start = std::chrono::system_clock::now();

	ofstream classout;
	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}

	string tree_type = "fist_tree";
	if (parameter.useHMT == 1) {
		tree_type = "hmt_tree";
	}

	idf = idf + "_" + tree_type;

	classout.open(CTOutputLocation + parameter.reachId + "_Prediction_selected_" + idf + ".txt");
	//classout.open(CTOutputLocation + parameter.fname + ".txt");

	for (int i = 0; i < mappredictions.size(); i++) {
		classout << mappredictions[i] << endl;
	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
		    if(mappredictions[index]>0)
		    {
		        prediction[row][col] = 1;
		    }
		    else
			    prediction[row][col] = mappredictions[index];
			index++;
		}
	}
	GDALDataset* srcDataset = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform[6];
	srcDataset->GetGeoTransform(geotransform);
	const OGRSpatialReference* poSRS = srcDataset->GetSpatialRef();
	;
	GeotiffWrite finalTiff((CTOutputLocation + parameter.reachId + "_Prediction_selected_" + idf + ".tif").c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform, poSRS);
	finalTiff.write(prediction);

	cout << "Selected Prediction finished!" << endl;
}

void cFlood::delta_prediction() {
	cout << "Delta prediction started!" << endl;

	auto start = std::chrono::system_clock::now();

	ofstream classout;
	classout.open(CTOutputLocation + parameter.reachId + "_Prediction_delta.txt");
	//classout.open(CTOutputLocation + parameter.fname + ".txt");

	for (int i = 0; i < mappredictions.size(); i++) {
		classout << mappredictions[i] << endl;
	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
			prediction[row][col] = mappredictions[index];
			index++;
		}
	}
	GDALDataset* srcDataset = (GDALDataset*)GDALOpen((CTInputLocation + CTFel).c_str(), GA_ReadOnly);
	double geotransform[6];
	srcDataset->GetGeoTransform(geotransform);
	const OGRSpatialReference* poSRS = srcDataset->GetSpatialRef();
	;
	GeotiffWrite finalTiff((CTOutputLocation + parameter.reachId + "_Prediction_delta.tif").c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform, poSRS);
	finalTiff.write(prediction);

	cout << "Delta Prediction finished!" << endl;
}

void cFlood::prediction() {
	cout << "prediction started!" << endl;


	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

		// no need to touch inferred regions
		if (data.rightInferredRegions[data.rightNodesInOrder[rightOrder]]) {
			continue;
		}

		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				if (data.inferedmaxCost_id2cost[data.reach_ids[rightOrder]] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					continue;
				}
				else {
					/*if (data.allNodes[nodid]->cost <= data.inferedmaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) {*/
					//if (data.allNodes[nodid]->fel <= data.inferedmaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) { // CHECK
					if (data.allNodes[nodid]->cost <= data.inferedmaxCost_id2cost[data.reach_ids[rightOrder]] && data.allNodes[nodid]->isNa == 0) { // CHECK NC
						// mappredictions[data.allNodes[nodid]->originalId] = 1;
						mappredictions[data.allNodes[nodid]->originalId] = data.rightNodesInOrder[rightOrder];
					}
					/*else if(data.allNodes[nodid]->isNa == 0){*/
					else if (mappredictions[data.allNodes[nodid]->originalId] == 0) {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}
					else {
						mappredictions[data.allNodes[nodid]->originalId] = -1; // NC
					}
				}
			}
		}
	}

	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		mappredictions[data.allNodes[data.rightNodesInOrder[i]]->originalId] = 1;
	}


	cout << "prediction finished!" << endl;
}





void cFlood::prediction_FIST() {
	cout << "prediction started!" << endl;
	//mappredictions.resize(parameter.orgPixelSize, -1);

	// CHECK: do not reset to -1, just skip inferred region and predict only on regions whose cost we found after interpolation
	// TODO: comment this for HMT, uncomment for FIST
	//std::fill(mappredictions.begin(), mappredictions.end(), -1);


	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		/*if (data.hasObservedPixelsRight[rightOrder] && data.rightbfsOrder[rightOrder].size()>= PIXELLIMT) {
			continue;
		}*/

		////// TODO: comment this for FIST, uncomment for HMT
		//if (std::find(data.rightInferredRegions.begin(), data.rightInferredRegions.end(), data.reach_ids[rightOrder]) != data.rightInferredRegions.end())
		//	continue;

		// no need to touch inferred regions
		if (data.rightInferredRegions[data.reach_ids[rightOrder]]) {
			continue;
		}

		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				if (data.inferedmaxCostRight[rightOrder] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					continue;
					/*if(data.allNodes[nodid]->isNa == 0)
						mappredictions[data.allNodes[nodid]->originalId] = 1;*/
				}
				else {

					if (data.allNodes[nodid]->cost <= data.inferedmaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) { // CHECK NC
						// mappredictions[data.allNodes[nodid]->originalId] = 1;
						mappredictions[data.allNodes[nodid]->originalId] = data.reach_ids[rightOrder];
					}
					else if (data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
						mappredictions[data.allNodes[nodid]->originalId] = data.reach_ids[rightOrder];
					}
					else {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}

				}
			}
		}
	}
	for (int i = 0; i < data.reach_ids.size(); i++) {
		//mappredictions[data.reach_ids_orig[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.reach_ids[i]]->originalId] = 1;
	}

	// nodes in river(high prob pixels from UNet should also be flooded)
	for (int i = 0; i < data.river_ids.size(); i++) {
		//mappredictions[data.river_ids_orig[i]] = 1;  //reach ids are the lowest point in the river
		mappredictions[data.allNodes[data.river_ids[i]]->originalId] = 1;
	}

	cout << "prediction finished!" << endl;
}



void cFlood::interpolate() {
	cout << "interpolation started!" << endl;
	//profile table before interpolation
	ofstream profiletable;
	// 	data.combinedCost.resize(data.reach_ids.size(), 0);
	// 	data.avgCost.resize(data.reach_ids.size(), 0);
		//find the regions with infered cost values (i.e values that are not -1)
	vector<int> stops;


	ofstream profiletable_right;


	profiletable_right.open(CTOutputLocation + "ProfileTables/" + CTTestCase + "_ProfileTable_preInterpolation_right_lambda_" + parameter.lambda_str.str() + ".csv");
	//profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.fname + "_ProfileTable_preInterpolation.csv");
	profiletable_right << "SourceId" << "," << "Right Infered Cost" << endl;
	for (int index = 0; index < data.rightNodesInOrder.size(); index++) {
		profiletable_right << data.rightNodesInOrder[index] << ","
			<< data.inferedmaxCostRight[index] << endl;
	}
	profiletable_right.close();

	//for right bank.
	int current = 0;
	while (current < data.rightNodesInOrder.size()) {
		if (data.inferedmaxCostRight[current] == -1 && current == 0) {

			//find the first reach node with non -1 max cost value
			int index = -1;
			for (int j = 1; j < data.rightNodesInOrder.size(); j++) {
				if (data.inferedmaxCostRight[j] != -1) {
					index = j;
					break;
				}
			}
			if (index == -1) {
				break;
			}
			double value = data.inferedmaxCostRight[index];
			for (int i = 0; i < index; i++) {
				data.inferedmaxCostRight[i] = value;
				data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = value;
			}
			current = index;


		}
		else if (data.inferedmaxCostRight[current] != -1) {
			//two cases
				//case 1: there are n points in between next reach that has cost value
				//case 2: there is no next point
			//find index of next reach node that has cost value
			int index = -1;
			int count = 0;
			double value = data.inferedmaxCostRight[current];
			for (int j = current + 1; j < data.rightNodesInOrder.size(); j++) {
				if (data.inferedmaxCostRight[j] != -1) {
					index = j;
					break;
				}
				count++;
			}
			if (index == -1) {// case 2
				for (int i = current + 1; i < data.rightNodesInOrder.size(); i++) {
					data.inferedmaxCostRight[i] = value;
					data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = value;
				}
				current = data.rightNodesInOrder.size();
				break;
			}
			else if (count == 0 && index == current + 1) {
				current = index;
			}
			else {
				double interval = (data.inferedmaxCostRight[index] - value) / count;
				for (int i = current + 1; i < index; i++) {
					data.inferedmaxCostRight[i] = data.inferedmaxCostRight[(i - 1)] + interval;
					data.inferedmaxCost_id2cost[data.rightNodesInOrder[i]] = data.inferedmaxCostRight[(i - 1)] + interval;
				}
				current = index;
			}
		}

	}
	cout << "interpolation finished!" << endl;

}



void cFlood::removeLink(vector<int>& v, int removeID) {
	v.erase(std::find(v.begin(), v.end(), removeID));
}

int cFlood::find(struct subset subsets[], int i)
{
	// find root and make root as parent of i (path compression)
	if (subsets[i].parent != i)
		subsets[i].parent = find(subsets, subsets[i].parent);

	return subsets[i].parent;
}

// A function that does union of two sets of x and y
// (uses union by rank)
void cFlood::Union(struct subset subsets[], int x, int y)
{
	int xroot = find(subsets, x);
	int yroot = find(subsets, y);

	// Attach smaller rank tree under root of high rank tree
	// (Union by Rank)
	if (subsets[xroot].rank < subsets[yroot].rank)
		subsets[xroot].parent = yroot;
	else if (subsets[xroot].rank > subsets[yroot].rank)
		subsets[yroot].parent = xroot;

	// If ranks are same, then make one as root and increment
	// its rank by one
	else
	{
		subsets[yroot].parent = xroot;
		subsets[xroot].rank++;
	}
}

// added by Saugat
// construct splitTree for each region on both left and right banks
void cFlood::splitTree() {


	// construct splitTree on right bank
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		int curIdx, neighborIndex;
		int row, column;

		vector<int> highestVertex(data.rightbfsOrder[i].size());
		subsets = (struct subset*)malloc(data.rightbfsOrder[i].size() * sizeof(struct subset));
		for (size_t j = 0; j < data.rightbfsOrder[i].size(); j++) {
			subsets[j].parent = j;
			subsets[j].rank = 0;
			highestVertex[j] = j;
		}

		for (size_t l = 0; l < data.rightbfsOrder[i].size(); l++) {
			curIdx = data.costIndexPairRight[i][l].second;
			row = curIdx / parameter.COLUMN;
			column = curIdx % parameter.COLUMN;
			highestVertex[curIdx] = curIdx;

			int curIdx_orig = data.rightIndex2PixelId[i][curIdx];

			double h1 = data.sortedCostIndexRight[i][curIdx];

			// check all 8 neighbors
			for (int j = max(0, row - 1); j <= min(parameter.ROW - 1, row + 1); j++) {
				for (int k = max(0, column - 1); k <= min(parameter.COLUMN - 1, column + 1); k++) {
					neighborIndex = j * parameter.COLUMN + k; //25

					if (data.NA[neighborIndex] == false && neighborIndex != curIdx) { // skip NA neighbor
						// skip neighbor from different region
						if (neighborIndex >= data.sortedCostIndexRight[i].size()) continue;

						int neighborIndex_orig = data.rightIndex2PixelId[i][neighborIndex];

						double h2 = data.sortedCostIndexRight[i][neighborIndex];

						if (h1 > h2) {
							int neighComponentID = find(subsets, neighborIndex);
							int currentComponetID = find(subsets, curIdx);
							if (neighComponentID == currentComponetID) {  //this means same as root2 == root1 but we don't need to find root2  //idea if they have same room they will point to same lowest vertex
								continue;
							}
							int currentHighestNodeIdx = highestVertex[neighComponentID];
							Union(subsets, curIdx, neighborIndex);

							int currentHighestNodeIdx_orig = data.rightIndex2PixelId[i][currentHighestNodeIdx];

						

							data.allNodes.at(currentHighestNodeIdx_orig)->childrenID.push_back(curIdx_orig);
							data.allNodes.at(curIdx_orig)->parentsID.push_back(currentHighestNodeIdx_orig); // TODO: check by reversing this


							/*data.allNodes.at(currentHighestNodeIdx_orig)->childrenID.push_back(curIdx_orig);
							data.allNodes.at(curIdx_orig)->parentsID.push_back(currentHighestNodeIdx_orig);*/

							int newComponentID = find(subsets, curIdx);
							highestVertex[newComponentID] = curIdx;
						}
					}
				}
			}
		}
	}



	// get new BFS order from split tree and then validate
	getNewBFSOrder();



}

void cFlood::getNewBFSOrder() {
	// clear the old bfs order
	data.leftbfsOrder.clear();
	data.rightbfsOrder.clear();


	//get root node for each tree
	// cout << "adj reach size: " << data.AdjustedReachNodes.size() << endl;

	// Left
	data.leftbfsRootNodes.resize(data.leftNodesInOrder.size(), -1); //leaf nodes with no children
	for (int i = 0; i < data.leftNodesInOrder.size(); i++) {
		int nodeId = data.leftNodesInOrder[i];
		//finding root in right tree
		int leftnode = -1;

		for (int j = 0; j < data.allNodes[nodeId]->childrenID.size(); j++) {
			int cid = data.allNodes[nodeId]->childrenID[j];
			leftnode = cid; //// get the first right node
		}

		if (leftnode != -1) {
			while (data.allNodes[leftnode]->childrenID.size() != 0) {
				leftnode = data.allNodes[leftnode]->childrenID[0]; //// get the right root node
			}
		}
		data.leftbfsRootNodes[i] = leftnode;

	}

	//get bfs order for each tree

	vector<int> bfsVisitedNew;
	bfsVisitedNew.resize(data.allNodes.size(), 0);

	for (int i = 0; i < data.leftNodesInOrder.size(); i++) {
		if (data.leftbfsRootNodes[i] == -1) {
			data.leftbfsOrder.push_back({});
		}
		else {
			data.leftbfsOrder.push_back(getBFSOrder(data.leftbfsRootNodes[i], bfsVisitedNew, 2));
		}
	}

	// Right
	data.rightbfsRootNodes.resize(data.rightNodesInOrder.size(), -1); //leaf nodes with no children
	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		int nodeId = data.rightNodesInOrder[i];
		//finding root in right tree
		int rightnode = -1;

		for (int j = 0; j < data.allNodes[nodeId]->childrenID.size(); j++) {
			int cid = data.allNodes[nodeId]->childrenID[j];
			rightnode = cid; //// get the first right node
		}

		if (rightnode != -1) {
			while (data.allNodes[rightnode]->childrenID.size() != 0) {
				rightnode = data.allNodes[rightnode]->childrenID[0]; //// get the right root node
			}
		}
		data.rightbfsRootNodes[i] = rightnode;

	}

	//get bfs order for each tree

	// vector<int> bfsVisitedNew;
	// bfsVisitedNew.resize(data.allNodes.size(), 0);

	for (int i = 0; i < data.rightNodesInOrder.size(); i++) {
		if (data.rightbfsRootNodes[i] == -1) {
			data.rightbfsOrder.push_back({});
		}
		else {
			data.rightbfsOrder.push_back(getBFSOrder(data.rightbfsRootNodes[i], bfsVisitedNew, 2));
		}
	}

}

void cFlood::displayTree(int TreeId) {

	cout << "Hidden Markov Tree Right Bank" << endl;
	cout << "Parent --- Current Node ---Child" << endl;
	cout << data.rightbfsOrder[TreeId].size() << endl;
	for (size_t k = 0; k < data.rightbfsOrder[TreeId].size(); k++) {
		cout << "k: " << k << endl;
		int i = data.rightbfsOrder[TreeId][k];

		if (data.allNodes[i]->parentsID.size() > 0) {
			cout << "<";
			for (int j = 0; j < data.allNodes[i]->parentsID.size(); j++) {
				if (j + 1 != data.allNodes[i]->parentsID.size())
					cout << data.allNodes[i]->parentsID[j] << ",";
				else
					cout << data.allNodes[i]->parentsID[j] << ">    <-----";
			}
		}
		else {
			cout << "<NUll>    <-----";
		}
		cout << i << "------>     ";
		if (data.allNodes[i]->childrenID.size() > 0) {
			cout << "<";
			for (int j = 0; j < data.allNodes[i]->childrenID.size(); j++) {
				if (j + 1 != data.allNodes[i]->childrenID.size())
					cout << data.allNodes[i]->childrenID[j] << ",";
				else
					cout << data.allNodes[i]->childrenID[j] << ">" << endl << endl;
			}
		}
		else {
			cout << "<NUll>" << endl << endl;
		}

	}
}



void cFlood::validateTreeLeft() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost

	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		//int root_id = data.leftbfsRootNodes[leftOrder];

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// traverse from each leaf nodes to bfs root
		for (int j = 0; j < leafnodes.size(); j++) {
			int nodeId = leafnodes[j];
			int parentId = nodeId;

			int bfsRoot = nodeId;

			int children_size = data.allNodes[nodeId]->childrenID.size();

			if (children_size == 0) {
				if (nodeId != bfsRoot) {
					cout << "Error: No children of pixel Id: " << nodeId << "---Tree ID: " << leftOrder << endl;
				}
			}

			if (children_size > 1) {
				if (std::find(data.AdjustedReachNodes.begin(), data.AdjustedReachNodes.end(), nodeId) == data.AdjustedReachNodes.end()) {
					// reach id can have multiple children on left and right bank
					cout << "Error: Multiple children of pixel Id " << nodeId << " ---Tree ID: " << leftOrder << endl;
				}
			}

			while (data.allNodes[nodeId]->childrenID.size() != 0) {
				int child_id = data.allNodes[nodeId]->childrenID[0];

				double child_cost = data.allNodes[child_id]->cost;
				double parent_cost = data.allNodes[nodeId]->cost;

				if (child_cost < parent_cost) {
					cout << "Error: Child's cost less than Parent's cost" << endl;
				}
				nodeId = child_id;
			}
		}
	}


}

void cFlood::validateTreeRight() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost

	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		//int root_id = data.leftbfsRootNodes[leftOrder];

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int nodeId = data.rightbfsOrder[rightOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// traverse from each leaf nodes to bfs root
		for (int j = 0; j < leafnodes.size(); j++) {
			int nodeId = leafnodes[j];
			int parentId = nodeId;

			int bfsRoot = nodeId;

			int children_size = data.allNodes[nodeId]->childrenID.size();

			if (children_size == 0) {
				if (nodeId != bfsRoot) {
					cout << "Error: No children of pixel Id: " << nodeId << "---Tree ID: " << rightOrder << endl;
				}
			}

			if (children_size > 1) {
				if (std::find(data.AdjustedReachNodes.begin(), data.AdjustedReachNodes.end(), nodeId) == data.AdjustedReachNodes.end()) {
					// reach id can have multiple children on left and right bank
					cout << "Error: Multiple children of pixel Id " << nodeId << " ---Tree ID: " << rightOrder << endl;
				}
			}

			while (data.allNodes[nodeId]->childrenID.size() != 0) {
				int child_id = data.allNodes[nodeId]->childrenID[0];

				double child_cost = data.allNodes[child_id]->cost;
				double parent_cost = data.allNodes[nodeId]->cost;

				if (child_cost < parent_cost) {
					cout << "Error: Child's cost less than Parent's cost" << endl;
				}
				nodeId = child_id;
			}
		}
	}
}

void cFlood::validateTreeInferenceLeftFIST() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost

	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		//int root_id = data.leftbfsRootNodes[leftOrder];

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->childrenID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// traverse from each leaf nodes to reach node
		for (int j = 0; j < leafnodes.size(); j++) {
			int nodeId = leafnodes[j];



			while (data.allNodes[nodeId]->parentsID.size() != 0) {
				int parent_id = data.allNodes[nodeId]->parentsID[0];

				/*double child_cost = data.allNodes[child_id]->cost;
				double parent_cost = data.allNodes[nodeId]->cost;

				if (child_cost < parent_cost) {
					cout << "Error: Child's cost less than Parent's cost" << endl;
				}*/

				// TODO: check infered flood here
				int child_class = mappredictions[data.allNodes[nodeId]->originalId];
				int parent_class = mappredictions[data.allNodes[parent_id]->originalId];

				if (child_class > parent_class) {
					cout << "Error on Left Bank: Child is flooded, parent is not" << " -- Reach Id: " << data.reach_ids[leftOrder] << endl;
				}

				nodeId = parent_id;
			}
		}
	}
}

void cFlood::validateTreeInferenceRightFIST() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost

	//for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
	//	//int root_id = data.leftbfsRootNodes[leftOrder];

	int rightOrder = 10453;

	vector<int> leafnodes;
	//step 1: get list of leaf nodes
	for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
		int nodeId = data.rightbfsOrder[rightOrder][i];
		if (data.allNodes[nodeId]->childrenID.size() == 0) {
			leafnodes.push_back(nodeId);

		}
	}

	map<int, bool> already_added;
	int sum_pixels = 0;

	// traverse from each leaf nodes to reach node
	for (int j = 0; j < leafnodes.size(); j++) {
		int nodeId = leafnodes[j];



		if (already_added.find(nodeId) == already_added.end()) {
			sum_pixels++;
			already_added.insert(make_pair(nodeId, true));
		}

		while (data.allNodes[nodeId]->parentsID.size() != 0) {


			int parent_id = data.allNodes[nodeId]->parentsID[0];

			if (already_added.find(parent_id) == already_added.end()) {
				sum_pixels++;
				already_added.insert(make_pair(parent_id, true));
			}



			nodeId = parent_id;
		}

		/*if (data.allNodes[nodeId]->originalId != data.reach_ids_orig[rightOrder]) {
			cout << " missing path from leaves to parent" << endl;
			break;
		}*/
	}

	cout << "sum all pixels: " << sum_pixels << endl;
	cout << "Region size: " << data.rightbfsOrder[rightOrder].size() << endl;
}
//}



void cFlood::validateTreeInferenceLeft() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost
	cout << "Validating tree inference left" << endl;
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		//int root_id = data.leftbfsRootNodes[leftOrder];

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// traverse from each leaf nodes to bfs root
		for (int j = 0; j < leafnodes.size(); j++) {
			int nodeId = leafnodes[j];
			int parentId = nodeId;

			int bfsRoot = nodeId;

			int children_size = data.allNodes[nodeId]->childrenID.size();

			if (children_size == 0) {
				if (nodeId != bfsRoot) {
					cout << "Error: No children of pixel Id: " << nodeId << "---Tree ID: " << leftOrder << endl;
				}
			}

			if (children_size > 1) {
				if (std::find(data.AdjustedReachNodes.begin(), data.AdjustedReachNodes.end(), nodeId) == data.AdjustedReachNodes.end()) {
					// reach id can have multiple children on left and right bank
					cout << "Error: Multiple children of pixel Id " << nodeId << " ---Tree ID: " << leftOrder << endl;
				}
			}

			while (data.allNodes[nodeId]->childrenID.size() != 0) {
				int child_id = data.allNodes[nodeId]->childrenID[0];

				double child_cost = data.allNodes[child_id]->cost;
				double parent_cost = data.allNodes[nodeId]->cost;

				if (child_cost < parent_cost) {
					cout << "Error: Child's cost less than Parent's cost" << endl;
				}

				// TODO: check infered flood here
				int child_class = mappredictions[data.allNodes[child_id]->originalId];
				int parent_class = mappredictions[data.allNodes[nodeId]->originalId];

				if (child_class > parent_class) {
					cout << "Error on Left Bank: Child is flooded, parent is not" << " -- Reach Id: " << data.reach_ids[leftOrder] << endl;
					//break;
				}

				nodeId = child_id;
			}
		}
	}
}

void cFlood::validateTreeInferenceRight() { // root_id means BFS root(from top to bottom) --> highest cost to lowest cost
	cout << "Validating tree inference right" << endl;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		//int root_id = data.leftbfsRootNodes[leftOrder];

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int nodeId = data.rightbfsOrder[rightOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// traverse from each leaf nodes to bfs root
		for (int j = 0; j < leafnodes.size(); j++) {
			int nodeId = leafnodes[j];
			int parentId = nodeId;

			int bfsRoot = nodeId;

			int children_size = data.allNodes[nodeId]->childrenID.size();

			if (children_size == 0) {
				if (nodeId != bfsRoot) {
					cout << "Error: No children of pixel Id: " << nodeId << "---Tree ID: " << rightOrder << endl;
				}
			}

			if (children_size > 1) {
				if (std::find(data.AdjustedReachNodes.begin(), data.AdjustedReachNodes.end(), nodeId) == data.AdjustedReachNodes.end()) {
					// reach id can have multiple children on left and right bank
					cout << "Error: Multiple children of pixel Id " << nodeId << " ---Tree ID: " << rightOrder << endl;
				}
			}

			while (data.allNodes[nodeId]->childrenID.size() != 0) {
				int child_id = data.allNodes[nodeId]->childrenID[0];

				double child_cost = data.allNodes[child_id]->cost;
				double parent_cost = data.allNodes[nodeId]->cost;

				if (child_cost < parent_cost) {
					cout << "Error: Child's cost less than Parent's cost" << endl;
				}

				// TODO: check infered flood here
				int child_class = mappredictions[data.allNodes[child_id]->originalId];
				int parent_class = mappredictions[data.allNodes[nodeId]->originalId];

				if (child_class > parent_class) {
					cout << "Error on Right Bank: Child is flooded, parent is not" << " -- Reach Id: " << data.reach_ids[rightOrder] << endl;

				}

				nodeId = child_id;
			}
			cout << "out of while right" << endl;
		}
		cout << "out of for right" << endl;
	}
	cout << "out of outer for right" << endl;
}

void cFlood::getStatistics() {
	ofstream stat_table;

	stat_table.open(CTOutputLocation + "Statistics\\" + parameter.reachId + "_Statistics.csv");
	/*stat_table << "Reach Id" << "," << "Pixel Id" << "," << "Class" << "," << "Bank" << "," << "No. of Children" << "," << "Multiple Children?" << "," << "Cost" << "," << "parent_id" << "," << "child_id" << endl;*/
	stat_table << "Reach Id" << "," << "Pixel Id" << "," << "Bank" << "," << "No. of Parents" << "," << "Multiple Parent?" << "," << "Cost" << "," << "child_id" << "," << "parent_id" << endl;

	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		int reach_id = data.reach_ids[i];
		bool multiple_parent;

		for (int j = 0; j < data.rightbfsOrder[i].size(); j++) {
			int pixelId = data.rightbfsOrder[i][j];
			int num_children = data.allNodes[pixelId]->childrenID.size();
			int num_parent = data.allNodes[pixelId]->parentsID.size();
			int bank = data.allNodes[pixelId]->bank;
			multiple_parent = false;
			if (num_parent > 1) multiple_parent = true;

			float cost = data.allNodes[pixelId]->cost;

			// get child
			int child_id = 0;
			if (data.allNodes[pixelId]->childrenID.size() > 0) child_id = data.allNodes[pixelId]->childrenID[0];

			// get children
			string parents_id = "";
			for (int k = 0; k < num_parent; k++) {
				int parent_id = data.allNodes[pixelId]->parentsID[k];
				parents_id = parents_id + "," + to_string(parent_id);
			}

			//int node_class = mappredictions[data.allNodes[pixelId]->originalId];

			// write to csv file
			//stat_table << reach_id << "," << pixelId << "," << node_class << "," << bank << "," << num_children << "," << multiple_children << "," << cost << "," << parent_id << "," << children_id << endl;
			stat_table << reach_id << "," << pixelId << "," << bank << "," << num_children << "," << multiple_parent << "," << cost << "," << child_id << "," << parents_id << endl;
		}
	}
	stat_table.close();
}



void cFlood::inference() {
	vector<int> inferVisited(parameter.allPixelSize, 0);




	//	//for right
	inferVisited.clear();
	inferVisited.resize(parameter.allPixelSize, 0);
	for (int i = 0; i < parameter.allPixelSize; i++) {
		data.allNodes[i]->correspondingNeighbour.clear();
		data.allNodes[i]->correspondingNeighbourClassOne.clear();
		data.allNodes[i]->correspondingNeighbourClassZero.clear();
	}
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		} // NC

		//if (data.rightbfsOrder[rightOrder].size() < PIXELLIMT) { // NC
		//	continue;
		//}



		int bfsTraversalOrderSize = (int)data.rightbfsOrder[rightOrder].size();
		//int bfsTraversalOrderSize = (int)data.bfsTraversalOrder.size();
		for (int node = bfsTraversalOrderSize - 1; node >= 0; node--) {
			int cur_node_id = data.rightbfsOrder[rightOrder][node];
			//int cur_node_id = data.bfsTraversalOrder[node];
			vector<int> rightbankchildrenID;
			for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				int child = data.allNodes[cur_node_id]->childrenID[c];
				/*if (data.allNodes[child]->bank == 2) {
					rightbankchildrenID.push_back(child);
				}*/

				rightbankchildrenID.push_back(child); // NC

			}
			//			data.allNodes[cur_node_id]->fi_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size()* cNum, 0);
			data.allNodes[cur_node_id]->fi_ChildList.resize(rightbankchildrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[cur_node_id]->fi[cls] = 0;
				data.allNodes[cur_node_id]->fo[cls] = 0;
			}

			//first figure out which neighbor fmessage passes to from current node pass n->? foNode;
			//idea: In bfs traversal order leave to root, check if next the node in bfs order is parent or child of the current node (should be child or parent of the current node)
			int foNode = -1;
			bool foNode_isChild = false;
			for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {  //check in parent list if found respective parent node is foNode
				int pid = data.allNodes[cur_node_id]->parentsID[p];
				if (!inferVisited[pid]) {
					foNode = pid;
					break;
				}
			}
			if (foNode == -1) {
				for (int c = 0; c < rightbankchildrenID.size(); c++) {
					int cid = rightbankchildrenID[c];
					if (!inferVisited[cid]) {
						foNode = cid;
						foNode_isChild = true;
						break;
					}
				}
			}
			data.allNodes[cur_node_id]->foNode = foNode;
			data.allNodes[cur_node_id]->foNode_ischild = foNode_isChild;
			if (cur_node_id == data.rightbfsRootNodes[rightOrder] && rightbankchildrenID.size() == 0) { //need to verify this changed && to || for IRONFIST project
				foNode_isChild = true;
			}

			//incoming message from visited child
			if (data.allNodes[cur_node_id]->childrenID.size() > 0) {

				for (int c = 0; c < rightbankchildrenID.size(); c++) {
					int child_id = rightbankchildrenID[c];

					if (child_id == foNode) {
						continue;
					}
					data.allNodes[cur_node_id]->correspondingNeighbour.push_back(child_id);
					for (int p = 0; p < data.allNodes[child_id]->parentsID.size(); p++) {
						int pid = data.allNodes[child_id]->parentsID[p];
						if (pid != cur_node_id) {
							data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
						}

					}
					vector<int> parentOfChildExcept_currentNode;
					for (int en = 0; en < data.allNodes[child_id]->parentsID.size(); en++) {
						if (data.allNodes[child_id]->parentsID[en] != cur_node_id) {
							parentOfChildExcept_currentNode.push_back(data.allNodes[child_id]->parentsID[en]);
						}

					}
					for (int cls = 0; cls < cNum; cls++) {  //cls represents current node class
						//double sumAccumulator = eln(0);   //should be 0 since we are summing it up//eln(1);//need to confirm
						double max = eln(0);
						vector<int> maxCorrespondingNeighbour;
						for (int c_cls = 0; c_cls < cNum; c_cls++) { //c_cls reperesnets child class label   Yc
							int max_bitCount = 1 << parentOfChildExcept_currentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent and child class label(given by c_cls)
								double productAccumulator = data.allNodes[child_id]->fo[c_cls];  //product with fo(c)
								vector<int>neighbourClass;
								neighbourClass.push_back(c_cls);
								int parentClsProd = 1; //p(c), product of parent classes for child c
								for (int p = 0; p < parentOfChildExcept_currentNode.size(); p++) {//calculating Product(fo(p)) for all parent of current child except the current node
									int pid = parentOfChildExcept_currentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									neighbourClass.push_back(parentClsValue);
									productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
								}
								//multiplying P(Yc|Ypc)
								parentClsProd *= cls;
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[c_cls][parentClsProd]);
								if (max < productAccumulator) {
									max = productAccumulator;
									maxCorrespondingNeighbour = neighbourClass;
								}
							}
						}
						data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)] = max;
						if (cls == 0) {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxCorrespondingNeighbour[t]);
							}
						}
						else {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxCorrespondingNeighbour[t]);
							}
						}
					}
				}
			}

			if (foNode_isChild) {  //means the current node has all visited parents
				if (data.allNodes[cur_node_id]->parentsID.size() == 0) {
					for (int cls = 0; cls < cNum; cls++) {
						data.allNodes[cur_node_id]->fi[cls] = parameter.elnPz[cls];
					}
				}
				else {
					for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
						int pid = data.allNodes[cur_node_id]->parentsID[p];
						data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
					}
					for (int cls = 0; cls < cNum; cls++) {
						double max = eln(0);
						vector<int> maxNeighbourClass;
						int max_bitCount = 1 << data.allNodes[cur_node_id]->parentsID.size();
						for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label
							vector<int> parentClass;
							double productAccumulator = eln(1);
							int parentClsProd = 1;
							for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
								int pid = data.allNodes[cur_node_id]->parentsID[p];
								int parentClsValue = (bitCount >> p) & 1;
								parentClass.push_back(parentClsValue);
								parentClsProd *= parentClsValue;
								productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
							}
							productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[cls][parentClsProd]);
							if (max < productAccumulator) {
								max = productAccumulator;
								maxNeighbourClass = parentClass;
							}
							//sumAccumulator = elnsum(sumAccumulator, productAccumulator);
						}
						data.allNodes[cur_node_id]->fi[cls] = max;
						if (cls == 0) {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxNeighbourClass[t]);
							}
						}
						else {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxNeighbourClass[t]);
							}
						}
					}
				}

				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < rightbankchildrenID.size(); c++) {
						int child_id = rightbankchildrenID[c];
						if (child_id == foNode) {
							continue;
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)]); //multiplying with al the child fi except the outgoing child
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi[cls]);  // multiplying with fi(n)_parent

					// TODO: CHECK THIS; we had elnPzn_xn previously d/t delta result
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);

					//productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id * cNum + cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}

			}

			else {  //message pass n-> parent there is no fi(n)_parent   //computes for root node as well
				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < rightbankchildrenID.size(); c++) {
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)]); //multiplying with al the child fi except the outgoing child
					}

					// TODO: CHECK THIS; we had elnPzn_xn previously d/t delta result
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);

					//productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[(cur_node_id*cNum + cls)]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}
			}

			inferVisited[cur_node_id] = 1;
		}
	}

	if (parameter.useHMT) {
		updateMapPrediction_right_hmt(); // for HMT tree(Split Tree)
	}
	else {
		updateMapPrediction_right(); // for FIST tree
		//updateMapPrediction_right_verify();
	}

	//verify_deltaResult_right();


	//cout << "validate Right Tree started" << endl;
	//validateTreeInferenceRight();
	//cout << "validate Right Tree ended" << endl;
}

//  function to calc std
float calc_std(vector<float> boundaryCostList, float sum, int n) {
	float avg = sum / n;

	float stdev = 0;
	float cost = 0;
	for (int i = 0; i < n; i++) {
		cost = boundaryCostList[i];
		stdev += pow(cost - avg, 2);
	}
	return sqrt(stdev / n);
}

// function to calc q3
float calc_q3(vector<float> boundaryCostList) {
	std::sort(boundaryCostList.begin(), boundaryCostList.end());
	int nn = boundaryCostList.size();
	//cout << "nn size: " << nn << endl;

	/*int q3_idx = round(3 * nn / 4);*/

	int q3_idx = round(nn / 2);
	//cout << "q3_idx: " << q3_idx << endl;
	//cout << "boundary: " << boundaryCostList[q3_idx] << endl;
	return boundaryCostList[q3_idx];
}

void cFlood::updateMapPrediction_left() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Leftbank: " << endl;

	// TODO: check
	ofstream boundary_costs_left;
	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}
	boundary_costs_left.open(CTOutputLocation + parameter.reachId + "_boundary_left_" + idf + ".txt");

	ofstream class_table;

	class_table.open(CTOutputLocation + "Statistics\\" + parameter.reachId + "_node_class.csv");
	class_table << "Reach_Id" << "," << "Pixel_Id" << "," << "Class_Original" << "," << "Class_Inferred" << endl;


	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}

		// Added by Saugat: filter out very far regions
		if (extra.leftRegions[leftOrder]->regionTooFar == true) {
			cout << "Region: " << data.reach_ids[leftOrder] << " too far!!!" << endl;
			continue;
		}

		// add reach ids of inferred regions
		data.leftInferredRegionsOld.push_back(data.reach_ids[leftOrder]);

		boundary_costs_left << "Reach Id: " << data.reach_ids[leftOrder] << " || "; // TODO: check

		ofstream boundaryTableLeft;

		vector<int> updatedNodes;

		string idf = "no_cutoff";
		if (parameter.useCutoff == 1) {
			idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
		}

		string tree_type = "fist_tree";
		if (parameter.useHMT == 1) {
			tree_type = "hmt_tree";
		}

		idf = idf + "_" + tree_type;

		int region_id = data.reach_ids[leftOrder];
		boundaryTableLeft.open(CTOutputLocation + "BoundaryTables\\" + "Left\\" + to_string(region_id) + "_" + idf + ".csv");
		boundaryTableLeft << "SourceId" << "," << "Cost" << endl;

		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[leftOrder] << "->" << data.leftbfsOrder[leftOrder].size() << "->"; // TODO: UNCOMMENT
		int bfsroot = data.leftbfsRootNodes[leftOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}

			//// TODO: added this to check result before inundation
			//if (nodeCls == 1) {
			//	mappredictions[data.allNodes[bfsroot]->originalId] = data.reach_ids[leftOrder];
			//}
			//else {
			//	mappredictions[data.allNodes[bfsroot]->originalId] = 0;
			//}




			//if (data.allNodes[bfsroot]->isPits == 1) {
			//	//cout << "Pits original";
			//	if (data.allNodes[bfsroot]->p < 0.5) {
			//		cout << "Pits original non flood";
			//		nodeCls = 0;
			//	}
			//	/*else if (data.allNodes[bfsroot]->p == 0.5) {
			//		cout << "Pits original not sure";
			//	}
			//	else {
			//		cout << "Pits original flood";
			//	}*/
			//}


			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;
			updatedNodes.push_back(data.allNodes[bfsroot]->originalId);

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}


					//if (data.allNodes[neigh_id]->isPits == 1) {
					//	//cout << "Pits original 2";
					//	if (data.allNodes[neigh_id]->p < 0.5) {
					//		cout << "Pits original non flood 2";
					//		cClass = 0;
					//	}
					//	/*else if (data.allNodes[neigh_id]->p == 0.5) {
					//		cout << "Pits original not sure 2";
					//	}
					//	else {
					//		cout << "Pits original flood 2";
					//	}*/
					//}

					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[leftOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;
					updatedNodes.push_back(data.allNodes[neigh_id]->originalId);

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}

				// TODO: check
				class_table << data.reach_ids[leftOrder] << "," << node << "," << data.allNodes[node]->isObserved << "," << nodeClass[node] << endl;

			}
		}

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->childrenID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
			/*if (data.allNodes[nodeId]->parentsID.size() > 1) {
				std::cout << data.allNodes[nodeId]->parentsID.size() << endl;
			}*/
		}

		//step 2: get chain length
		for (int i = 0; i < leafnodes.size(); i++) {
			int chainLength = 0;
			int nodeId = leafnodes[i];
			int leafnodeId = nodeId;
			int roughness = 1;
			float boundary_cost = -1.0;
			float chainMaxCost = data.allNodes[nodeId]->cost;
			pair<float, float> temp_maxCost_boundaryCost_pair/* = make_pair(chainMaxCost, -1.0)*/;
			pair<int, float> temp_cl_cost_pair /*= make_pair(0, -1.0)*/;
			pair<int, int> temp_chainLength_id/* = make_pair(-1, -1)*/;
			pair<float, int> temp_chainMaxCost_id/* = make_pair(chainMaxCost, -1)*/;
			pair<float, int> temp_cost_boundaryNode_pair /*= make_pair(-1.0,-1)*/;
			pair<int, vector<int>> temp_boundaryNode_leafNodes_pair;
			if (EB == 0) {
				temp_cl_cost_pair = make_pair(0, -1.0);
				temp_chainLength_id = make_pair(-1, -1);
			}
			else if (EB == 1) {
				temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);
				temp_chainMaxCost_id = make_pair(chainMaxCost, -1);
			}
			else if (EB == 2) {
				temp_cost_boundaryNode_pair = make_pair(-1.0, -1);
			}
			//pair<float, float> temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);

			while (data.allNodes[nodeId]->parentsID.size() != 0) {
				if (nodeClass[nodeId] == 1 && boundary_cost < 0.0) { // finding first flood boundary node
					// TODO: remove this
					//if (data.reach_ids[leftOrder] == 3050550) cout << "3050550 first flood boundary node found" << endl;

					bool allchildObserved = true;
					bool observed = true;
					if (BOUNDARY_NODES_OBSERVED == 1) {  // 0 does not exclude any branches// 1 consider pits and tree and unobserved and exclude those branches
						//2 excludes branches if the boundary of flood and dry is overlapping with pits layer
						for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
							int childid = data.allNodes[nodeId]->childrenID[idx];
							if (data.allNodes[childid]->isObserved != 1) {
								allchildObserved = false;
								break;
							}
						}

						if (data.allNodes[nodeId]->isObserved == 0) {
							observed = false;
							break;
						}
					}
					if (BOUNDARY_NODES_OBSERVED == 2) {  // uncomment after adding pits and Na identifiers
						for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
							int childid = data.allNodes[nodeId]->childrenID[idx];
							if (data.allNodes[childid]->isPits == 1 || data.allNodes[nodeId]->isNa == 1) {
								allchildObserved = false;
								break;
							}
						}

						if (data.allNodes[nodeId]->isNa == 1 || data.allNodes[nodeId]->isPits == 1) {
							/*if (data.reach_ids[leftOrder] == 3618868)
								cout << "observed false for reach : " << data.reach_ids[leftOrder];*/
							observed = false;
						}
					}
					if (data.allNodes[nodeId]->roughness == 1 || !observed || !allchildObserved || data.allNodes[nodeId]->isNa == 1) {
						/*if(data.reach_ids[leftOrder] == 3618868)
							cout << "breaking for reach: " << data.reach_ids[leftOrder];*/
						break;
					}
					boundary_cost = data.allNodes[nodeId]->cost;
					if (EB == 0) {
						temp_cl_cost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainLength_id.second = leafnodeId;
					}
					else if (EB == 1) {
						temp_maxCost_boundaryCost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainMaxCost_id.second = leafnodeId;
					}
					else if (EB == 2) {
						temp_cost_boundaryNode_pair.first = data.allNodes[nodeId]->cost;
						temp_cost_boundaryNode_pair.second = nodeId;

						if (data.boundaryNode_leafNodes_Map.find(nodeId) == data.boundaryNode_leafNodes_Map.end()) {
							vector<int> leafs;
							leafs.push_back(leafnodeId);
							data.boundaryNode_leafNodes_Map.insert(make_pair(nodeId, leafs));
						}
						else {
							data.boundaryNode_leafNodes_Map[nodeId].push_back(leafnodeId);
						}
					}
					data.leafNode_boundaryNodes.insert(make_pair(leafnodeId, nodeId));
					roughness = 0;
					if ((EB == 1) || (EB == 2)) {
						break;
					}
				}
				if (EB == 0) {
					chainLength++;
				}
				nodeId = data.allNodes[nodeId]->parentsID[0];
			}
			if (EB == 0) {
				temp_cl_cost_pair.first = chainLength;
				temp_chainLength_id.first = chainLength;
			}

			// TODO: remove
			float b_c = temp_cost_boundaryNode_pair.first;
			//if (data.reach_ids[leftOrder] == 1046162) cout << "1046162 boundary cost: " << b_c << endl;

			//if (data.reach_ids[leftOrder] == 3050550) cout << "Reach Id: " << data.reach_ids[leftOrder] << " Roughness: " << roughness << endl;

			if (roughness == 0) {
				//if (data.reach_ids[leftOrder] == 3050550) cout << "3050550 roughness 0" << endl;
				if (EB == 0) {
					data.chainLength_cost_Pairs.push_back(temp_cl_cost_pair);
					data.chainLength_nodeid_Pairs.push_back(temp_chainLength_id);
				}
				else if (EB == 1) {
					data.maxChainCost_cost_Pairs.push_back(temp_maxCost_boundaryCost_pair);
					data.maxCost_nodeid_Pairs.push_back(temp_chainMaxCost_id);
				}
				else if (EB == 2) {
					data.cost_boundaryNode_Pairs.insert(temp_cost_boundaryNode_pair);
				}
			}
		}
		if (EB == 0) {
			if (data.chainLength_cost_Pairs.size() != data.chainLength_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}
		else if (EB == 1) {
			if (data.maxChainCost_cost_Pairs.size() != data.maxCost_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}

		if (EB == 0) {
			// using chain length or max cost
			if (data.chainLength_cost_Pairs.size() != 0) {
				sort(data.chainLength_cost_Pairs.rbegin(), data.chainLength_cost_Pairs.rend());
				sort(data.chainLength_nodeid_Pairs.rbegin(), data.chainLength_nodeid_Pairs.rend());

				//top 20 percent
				//int top = data.chainLength_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.chainLength_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.chainLength_cost_Pairs[j].second;
					if (data.chainLength_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for chain length
					data.boundary_LeafNodes.push_back(data.chainLength_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostLeft[leftOrder] = avg;
				}

				std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
			}
		}

		else if (EB == 1) {
			// using max cost branches
			if (data.maxChainCost_cost_Pairs.size() != 0) {

				sort(data.maxChainCost_cost_Pairs.rbegin(), data.maxChainCost_cost_Pairs.rend());
				sort(data.maxCost_nodeid_Pairs.rbegin(), data.maxCost_nodeid_Pairs.rend());


				//top 20 percent
				//int top = data.maxChainCost_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.maxChainCost_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.maxChainCost_cost_Pairs[j].second;
					if (data.maxCost_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for max cost
					data.boundary_LeafNodes.push_back(data.maxCost_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostLeft[leftOrder] = avg;
				}
				if (EB == 1) {
					std::cout << data.maxChainCost_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
				}
				else {
					std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
				}

			}
		}

		else if (EB == 2) {
			//using boundary cost

			// maintain vector of costs
			vector<float> boundaryCostList;

			if (data.cost_boundaryNode_Pairs.size() != 0) {
				if (DEBUG_OUTPUT == 1) {
					vector<float> infered_cost;
					set<pair<float, int>>::reverse_iterator it;
					for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
						pair<float, int> p = *it;
						infered_cost.push_back(p.first);
					}
					//ofstream inferedCosts;
					//inferedCosts.open(CTOutputLocation + "ProfileTables\\"+ to_string(leftOrder) + "_"+ to_string(data.reach_ids[leftOrder])+"_"+parameter.reachId + "_" + "inferedCosts_left.txt");
					////classout.open(CTOutputLocation + CTPrediction);
					//for (int i = 0; i < infered_cost.size(); i++) {
					//	inferedCosts << infered_cost[i] << endl;
					//}
					//inferedCosts.close();
				}

				// TODO: don't use cutoff for experiment
				//int top = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
				//int top = data.cost_boundaryNode_Pairs.size();

				int top = data.cost_boundaryNode_Pairs.size();
				int top_20 = 0;
				int top_80 = 0;
				if (parameter.useCutoff == 1) {
					top_20 = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
					top_80 = data.cost_boundaryNode_Pairs.size() * (1 - parameter.cutoff_percentage);
					top = top_80 - top_20;
				}


				float sum = 0.0;
				set<pair<float, int>>::reverse_iterator it;
				int counter = 0;
				for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
					// throw away top 20 and bottom 20 percent
					if (parameter.useCutoff == 1) {
						counter++;
						if (counter <= top_20) continue;
						if (counter > top_80) break;
					}


					pair<float, int> p = *it;
					int bNode = p.second;
					sum = sum + p.first;

					boundaryCostList.push_back(p.first);

					// get children
					string children_id = "";
					for (int k = 0; k < data.allNodes[bNode]->childrenID.size(); k++) {
						int child_id = data.allNodes[bNode]->childrenID[k];
						children_id = children_id + "," + to_string(child_id);
					}

					// dump boundary cost to text file
					boundary_costs_left << bNode << "," << p.first << "#" << data.allNodes[bNode]->parentsID[0] << "!" << children_id << "@";

					// dump boundary cost of each region to separate csv file
					boundaryTableLeft << bNode << "," << p.first << endl;

					for (int n = 0; n < data.boundaryNode_leafNodes_Map[bNode].size(); n++) {
						int lNode = data.boundaryNode_leafNodes_Map[bNode][n];
						data.boundary_LeafNodes.push_back(lNode);
					}

					// don't use cutoff
					if (parameter.useCutoff != 1) {
						counter++;
						if (counter == top) {
							break;
						}
					}

				}
				float avg = -1.0;
				if (top != 0) {
					avg = sum / top;
				}

				// TODO: get std
				if (top != 0) {
					float stdev = calc_std(boundaryCostList, sum, top);

					cout << "Reach ID: " << data.reach_ids[leftOrder] << " std: " << stdev << endl;

					if (data.reach_ids[leftOrder] == 3050550) cout << "3050550 avg cost: " << avg << endl;

					if (avg > parameter.minCost) {
						data.inferedmaxCostLeft[leftOrder] = avg;
					}

					extra.standardDeviationLeft[leftOrder] = stdev;

					if (stdev > 2.0) {
						data.inferedmaxCostLeft[leftOrder] = -1;
						for (int mm = 0; mm < updatedNodes.size(); mm++) {
							mappredictions[updatedNodes[mm]] = -1; // reset updated predn to -1 // TODO: check
						}
					}
				}

				std::cout << "->" << data.inferedmaxCostLeft[leftOrder] << endl;

			}
		}

		data.chainLength_cost_Pairs.clear();
		data.chainLength_cost_Pairs.shrink_to_fit();
		data.chainLength_nodeid_Pairs.clear();
		data.chainLength_nodeid_Pairs.shrink_to_fit();

		data.maxChainCost_cost_Pairs.clear();
		data.maxChainCost_cost_Pairs.shrink_to_fit();
		data.maxCost_nodeid_Pairs.clear();
		data.maxCost_nodeid_Pairs.shrink_to_fit();

		data.cost_boundaryNode_Pairs.clear();
		data.boundaryNode_leafNodes_Map.clear();

		boundary_costs_left << endl;
		boundaryTableLeft.close();

	}
	boundary_costs_left.close();

	class_table.close();
	cout << "Leftbank inference finished" << endl;
}

void cFlood::updateMapPrediction_left_new() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Leftbank: " << endl;

	// TODO: check
	ofstream costs_left;
	costs_left.open(CTOutputLocation + parameter.reachId + "_costs_left.txt");

	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {

		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}

		// add reach ids of inferred regions
		data.leftInferredRegionsOld.push_back(data.reach_ids[leftOrder]);

		costs_left << "Reach Id: " << data.reach_ids[leftOrder] << " || "; // TODO: check

		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[leftOrder] << "->" << data.leftbfsOrder[leftOrder].size() << "->";
		int bfsroot = data.leftbfsRootNodes[leftOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) { // TODO: check eln prob
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}

			if (nodeCls == 1) {
				mappredictions[data.allNodes[bfsroot]->originalId] = data.reach_ids[leftOrder];
			}
			else {
				mappredictions[data.allNodes[bfsroot]->originalId] = 0;
			}

			//mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}
					nodeClass[neigh_id] = cClass;

					if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[leftOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}
					//mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}

		SUM_COST = 0.0;
		COUNTER = 0;
		double boundary_cost = -1.0;


		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// maintain a vector of nodeId whose cost we already summed up
		vector<int> visited_nodes;

		// first traverse from bfs root towards parents and get the flood boundary
		int root_id = bfsroot;
		while (data.allNodes[root_id]->parentsID.size() != 0) {
			if (nodeClass[root_id] == 1) {
				SUM_COST += data.allNodes[root_id]->cost;
				COUNTER++;

				visited_nodes.push_back(root_id);

				// get parents
				string parents_id = "";
				for (int k = 0; k < data.allNodes[root_id]->parentsID.size(); k++) {
					int p_id = data.allNodes[root_id]->parentsID[k];
					parents_id = parents_id + "," + to_string(p_id);
				}

				int children_id = -1;
				if (data.allNodes[root_id]->childrenID.size() != 0) {
					children_id = data.allNodes[root_id]->childrenID[0];
				}
				costs_left << root_id << "," << data.allNodes[root_id]->cost << "#" << parents_id << "!" << children_id << "@";

				// break the loop as soon as we find the first 1(flood node)
				break;
			}
			root_id = data.allNodes[root_id]->parentsID[0];
		}

		// traverse from every leaf nodes(lower branches) towards children to find flood boundary
		// search for first 0(dry) and consider its parent as flood boundary
		// if we find a 0, all above that branch should be 0
		// basically, find the last 1(flood node) on each branch from bottom to top and consider it as flood boundary
		for (int i = 0; i < leafnodes.size(); i++) {
			int nodeId = leafnodes[i];
			int parentId = nodeId;

			while (data.allNodes[nodeId]->childrenID.size() != 0) {
				if (nodeClass[nodeId] == 0) {
					if (nodeId != parentId) { // leaf node is dry so all above in this branch should be dry
						// skip nodes whose value we already summed up
						if (std::find(visited_nodes.begin(), visited_nodes.end(), parentId) != visited_nodes.end())
							break;

						SUM_COST += data.allNodes[parentId]->cost;
						COUNTER++;

						visited_nodes.push_back(parentId);

						// get parents
						string parents_id = "";
						for (int k = 0; k < data.allNodes[parentId]->parentsID.size(); k++) {
							int p_id = data.allNodes[parentId]->parentsID[k];
							parents_id = parents_id + "," + to_string(p_id);
						}

						costs_left << parentId << "," << data.allNodes[parentId]->cost << "#" << parents_id << "!" << nodeId << "@";
					}
					break;

				}
				parentId = nodeId;
				nodeId = data.allNodes[nodeId]->childrenID[0];
			}
		}

		if (COUNTER > 0) boundary_cost = SUM_COST / COUNTER;

		cout << "boundary cost: " << boundary_cost << endl;

		data.inferedmaxCostLeft[leftOrder] = boundary_cost;

		std::cout << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
		costs_left << endl;

	}
	costs_left.close();
	cout << "Leftbank inference finished" << endl;
}

void cFlood::updateMapPrediction_left_hmt() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Leftbank: " << endl;

	// TODO: check
	ofstream costs_left;
	costs_left.open(CTOutputLocation + parameter.reachId + "_costs_left.txt");

	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {

		//// TODO: remove this
		//if (data.reach_ids[leftOrder] != 30008661) {
		//	continue;
		//}

		float max_cost = 0;

		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}

		// add reach ids of inferred regions
		data.leftInferredRegions.insert(make_pair(data.leftNodesInOrder[leftOrder], true));

		ofstream boundaryTableLeft;

		string idf = "no_cutoff";
		if (parameter.useCutoff == 1) {
			idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
		}

		string tree_type = "fist_tree";

		if (parameter.useHMT == 1) {
			tree_type = "hmt_tree";
		}

		idf = idf + "_" + tree_type;

		int region_id = data.leftNodesInOrder[leftOrder];
		boundaryTableLeft.open(CTOutputLocation + "BoundaryTables/" + to_string(region_id) + "_Left" + idf + ".csv");
		boundaryTableLeft << "SourceId" << "," << "Cost" << endl;

		costs_left << "Reach Id: " << data.leftNodesInOrder[leftOrder] << " || "; // TODO: check

		float maxCost = MAXCOST;
		std::cout << endl << data.leftNodesInOrder[leftOrder] << "->" << data.leftbfsOrder[leftOrder].size() << "->";
		int bfsroot = data.leftbfsRootNodes[leftOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}

			/*if (nodeCls == 1) {
				mappredictions[data.allNodes[bfsroot]->originalId] = data.reach_ids[leftOrder];
			}
			else {
				mappredictions[data.allNodes[bfsroot]->originalId] = 0;
			}*/

			// CHECK: set NA regions as -1
			int nodeCls_new = nodeCls;
			/*if (data.allNodes[bfsroot]->isNa == 1) {
				nodeCls_new = -1;
			}*/

			/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls_new; // NC

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}
					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[leftOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					// CHECK: set NA regions as -1
					int nodeCls_new = nodeCls;
					/*if (data.allNodes[neigh_id]->isNa == 1) {
						nodeCls_new = -1;
					}*/

					/*mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;*/
					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls_new; // NC

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}


		SUM_COST = 0.0;
		COUNTER = 0;

		double boundary_cost = -1.0;
		double frontierNode = -1;

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// 		cout << "leafnodes collected: " << leafnodes.size() << endl;

				// maintain a vector of nodeId whose cost we already summed up
		vector<float> visited_nodes;
		map<int, bool> visited_nodes_map;




// 		cout << "traverse from reach node to root" << endl;
		int nodeId = data.leftNodesInOrder[leftOrder];
		int parentId = nodeId;
		double maxCostFrontierNode = -1;

		while (data.allNodes[nodeId]->childrenID.size() != 0) {
			/*if (data.allNodes[nodeId]->cost > max_cost && data.allNodes[nodeId]->isNa == 0) {
				max_cost = data.allNodes[nodeId]->cost;
			}*/

			if (data.allNodes[nodeId]->cost > max_cost && data.allNodes[nodeId]->isNa == 0 && nodeClass[nodeId] == 1) {
				max_cost = data.allNodes[nodeId]->cost;
				maxCostFrontierNode = nodeId;
			}

			if (nodeClass[nodeId] == 0 && data.allNodes[nodeId]->isNa == 0) {
				if (data.leftNodesInOrder[leftOrder] == 30008661) {
					cout << "nodeclass 0" << endl;
				}
				if (nodeId != parentId) { // leaf node is dry so all above in this branch should be dry
					// skip nodes whose value we already summed up
					/*if (std::find(visited_nodes.begin(), visited_nodes.end(), parentId) != visited_nodes.end())
						break*/;



						SUM_COST += data.allNodes[parentId]->cost;
						//SUM_COST += data.allNodes[parentId]->fel; // NC
						COUNTER++;
						frontierNode = parentId;



						boundaryTableLeft << parentId << "," << data.allNodes[parentId]->cost << endl;
				}
				break;

			}
			//else if (nodeClass[nodeId] == 1 && data.allNodes[nodeId]->isNa == 1) { // if we find a isNa node, all below that is flood for else case
			//	if (nodeId != parentId) {
			//		SUM_COST += data.allNodes[parentId]->cost;
			//		//SUM_COST += data.allNodes[parentId]->fel; // NC
			//		COUNTER++;

			//		boundaryTableLeft << parentId << "," << data.allNodes[parentId]->cost << endl;
			//	}
			//	break;
			//}
			parentId = nodeId;
			nodeId = data.allNodes[nodeId]->childrenID[0];
		}

		if (COUNTER > 0) boundary_cost = SUM_COST / COUNTER;

		// 		cout << "counter: " << COUNTER << endl;
		cout << "boundary cost: " << boundary_cost << endl;

		// if all the nodes are flooded for a region; get max cost among flooded
		if (boundary_cost == -1) {
			boundary_cost = max_cost;
			frontierNode = maxCostFrontierNode;
		}

		data.inferedmaxCostLeft[leftOrder] = boundary_cost;
		// 		data.inferredFloodFrontier[leftOrder] = frontierNode; 
		data.inferredFloodFrontier_regionId2Cost.insert(make_pair(data.leftNodesInOrder[leftOrder], boundary_cost));
		data.inferredFloodFrontier_regionId2nodeId.insert(make_pair(data.leftNodesInOrder[leftOrder], frontierNode));
		data.inferredFloodFrontier_regionId2Size.insert(make_pair(data.leftNodesInOrder[leftOrder], data.leftbfsOrder[leftOrder].size()));

		std::cout << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
		costs_left << endl;
		boundaryTableLeft.close();

	}
	costs_left.close();
	cout << "Leftbank inference finished" << endl;
	cout << "#################################" << endl;
}


void cFlood::updateMapPrediction_right_hmt() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Rightbank: " << endl;

	// TODO: check
	ofstream costs_right;
	costs_right.open(CTOutputLocation + parameter.reachId + "_costs_right.txt");

	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

		//// TODO: remove this
		//if (data.rightNodesInOrder[rightOrder] != 30008661) {
		//	continue;
		//}

		float max_cost = 0;

		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}

		// add reach ids of inferred regions
		data.rightInferredRegions.insert(make_pair(data.rightNodesInOrder[rightOrder], true));

		ofstream boundaryTableRight;

		string idf = "no_cutoff";
		if (parameter.useCutoff == 1) {
			idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
		}

		string tree_type = "fist_tree";

		if (parameter.useHMT == 1) {
			tree_type = "hmt_tree";
		}

		idf = idf + "_" + tree_type;

		int region_id = data.rightNodesInOrder[rightOrder];
		boundaryTableRight.open(CTOutputLocation + "BoundaryTables/" + to_string(region_id) + "_" + idf + ".csv");
		boundaryTableRight << "SourceId" << "," << "Cost" << endl;

		costs_right << "Reach Id: " << data.rightNodesInOrder[rightOrder] << " || "; // TODO: check

		float maxCost = MAXCOST;
		std::cout << endl << data.rightNodesInOrder[rightOrder] << "->" << data.rightbfsOrder[rightOrder].size() << "->";
		int bfsroot = data.rightbfsRootNodes[rightOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}

			/*if (nodeCls == 1) {
				mappredictions[data.allNodes[bfsroot]->originalId] = data.rightNodesInOrder[rightOrder];
			}
			else {
				mappredictions[data.allNodes[bfsroot]->originalId] = 0;
			}*/

			// CHECK: set NA regions as -1
			int nodeCls_new = nodeCls;
			/*if (data.allNodes[bfsroot]->isNa == 1) {
				nodeCls_new = -1;
			}*/

			/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls_new; // NC

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}
					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.rightNodesInOrder[rightOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					// CHECK: set NA regions as -1
					int nodeCls_new = nodeCls;
					/*if (data.allNodes[neigh_id]->isNa == 1) {
						nodeCls_new = -1;
					}*/

					/*mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;*/
					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls_new; // NC

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}


		SUM_COST = 0.0;
		COUNTER = 0;

		double boundary_cost = -1.0;
		double frontierNode = -1;

		vector<int> leafnodes;
		//step 1: get list of leaf nodes
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int nodeId = data.rightbfsOrder[rightOrder][i];
			if (data.allNodes[nodeId]->parentsID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}

		// 		cout << "leafnodes collected: " << leafnodes.size() << endl;

				// maintain a vector of nodeId whose cost we already summed up
		vector<float> visited_nodes;
		map<int, bool> visited_nodes_map;


// 		cout << "traverse from reach node to root" << endl;
		int nodeId = data.rightNodesInOrder[rightOrder];
		int parentId = nodeId;
		double maxCostFrontierNode = -1;

		while (data.allNodes[nodeId]->childrenID.size() != 0) {
			/*if (data.allNodes[nodeId]->cost > max_cost && data.allNodes[nodeId]->isNa == 0) {
				max_cost = data.allNodes[nodeId]->cost;
			}*/

			if (data.allNodes[nodeId]->cost > max_cost && data.allNodes[nodeId]->isNa == 0 && nodeClass[nodeId] == 1) {
				max_cost = data.allNodes[nodeId]->cost;
				maxCostFrontierNode = nodeId;
			}

			if (nodeClass[nodeId] == 0 && data.allNodes[nodeId]->isNa == 0) {
				if (data.rightNodesInOrder[rightOrder] == 30008661) {
					cout << "nodeclass 0" << endl;
				}
				if (nodeId != parentId) { // leaf node is dry so all above in this branch should be dry
					// skip nodes whose value we already summed up
					/*if (std::find(visited_nodes.begin(), visited_nodes.end(), parentId) != visited_nodes.end())
						break*/;

						//if (visited_nodes_map.find(parentId) == visited_nodes_map.end())
						/*if (visited_nodes_map[parentId])
							break;*/

						SUM_COST += data.allNodes[parentId]->cost;
						//SUM_COST += data.allNodes[parentId]->fel; // NC
						COUNTER++;
						frontierNode = parentId;



						boundaryTableRight << parentId << "," << data.allNodes[parentId]->cost << endl;
				}
				break;

			}

			parentId = nodeId;
			nodeId = data.allNodes[nodeId]->childrenID[0];
		}

		if (COUNTER > 0) boundary_cost = SUM_COST / COUNTER;

		// 		cout << "counter: " << COUNTER << endl;
		cout << "boundary cost: " << boundary_cost << endl;

		// if all the nodes are flooded for a region; get max cost among flooded
		if (boundary_cost == -1) {
			boundary_cost = max_cost;
			frontierNode = maxCostFrontierNode;
		}

		data.inferedmaxCostRight[rightOrder] = boundary_cost;
		// 		data.inferredFloodFrontier[rightOrder] = frontierNode; 
		data.inferredFloodFrontier_regionId2Cost.insert(make_pair(data.rightNodesInOrder[rightOrder], boundary_cost));
		data.inferredFloodFrontier_regionId2nodeId.insert(make_pair(data.rightNodesInOrder[rightOrder], frontierNode));
		data.inferredFloodFrontier_regionId2Size.insert(make_pair(data.rightNodesInOrder[rightOrder], data.rightbfsOrder[rightOrder].size()));

		std::cout << "->" << data.inferedmaxCostRight[rightOrder] << endl;
		costs_right << endl;
		boundaryTableRight.close();

	}
	costs_right.close();
	cout << "Rightbank inference finished" << endl;
}


void cFlood::updateMapPrediction_right_verify() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for right trees
	std::cout << "Rightbank: " << endl;

	// TODO: check
	//ofstream boundary_costs_right;
	//boundary_costs_right.open(CTOutputLocation + parameter.reachId + "_boundary_right.txt");

	//for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

	int rightOrder = 10453;


	ofstream boundaryTableRight;

	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}

	string tree_type = "fist_tree";

	if (parameter.useHMT == 1) {
		tree_type = "hmt_tree";
	}

	idf = idf + "_" + tree_type;

	int region_id = data.reach_ids[rightOrder];
	boundaryTableRight.open(CTOutputLocation + "BoundaryTables/" + to_string(region_id) + "_" + idf + ".csv");
	boundaryTableRight << "SourceId" << "," << "Cost" << endl;

	float maxCost = MAXCOST;
	//std::cout << endl << data.reach_ids[rightOrder] << "->" << data.rightbfsOrder[rightOrder].size() << "->"; // TODO: UNCOMMENT
	int bfsroot = data.rightbfsRootNodes[rightOrder];
	int orgId;
	if (bfsroot != -1) {
		queue<int> que;
		int nodeCls;
		if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
			nodeCls = 0;
		}
		else {
			nodeCls = 1;
		}



		// CHECK: set NA regions as -1
		int nodeCls_new = nodeCls;
		/*if (data.allNodes[bfsroot]->isNa == 1) {
			nodeCls_new = -1;
		}*/

		/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
		mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls_new; // NC

		//updatedNodes.push_back(data.allNodes[bfsroot]->originalId);

		nodeClass[bfsroot] = nodeCls;
		que.push(bfsroot);
		while (!que.empty()) {
			int node = que.front();
			que.pop();
			int nodeCls = nodeClass[node];
			visited[node] = 1;
			for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
				int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
				int cClass;
				if (nodeCls == 0) {
					cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
				}
				else {
					cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
				}

				nodeClass[neigh_id] = cClass;

				// CHECK: set NA regions as -1
				int nodeCls_new_neigh = nodeCls;
				/*if (data.allNodes[neigh_id]->isNa == 1) {
					nodeCls_new_neigh = -1;
				}*/

				/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
				mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls_new_neigh; // NC

				//updatedNodes.push_back(data.allNodes[neigh_id]->originalId);

				if (!visited[neigh_id]) {
					que.push(neigh_id);
				}

			}
		}
	}

	validateTreeInferenceRightFIST();
	//}
}


void cFlood::updateMapPrediction_right() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for right trees
	std::cout << "Rightbank: " << endl;

	// TODO: check
	//ofstream boundary_costs_right;
	//boundary_costs_right.open(CTOutputLocation + parameter.reachId + "_boundary_right.txt");

	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

		// TODO: CHECK NC
		////vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		//
		//
		//if (data.rightbfsOrder[rightOrder].size() < PIXELLIMT) { // NC
		//	continue;
		//}

	

		vector<int> updatedNodes;

		// add reach ids of inferred regions
		data.rightInferredRegions.insert(make_pair(data.reach_ids[rightOrder], true));

		//boundary_costs_right << "Reach Id: " << data.reach_ids[rightOrder] << " || "; // TODO: check

		ofstream boundaryTableRight;

		string idf = "no_cutoff";
		if (parameter.useCutoff == 1) {
			idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
		}

		string tree_type = "fist_tree";

		if (parameter.useHMT == 1) {
			tree_type = "hmt_tree";
		}

		idf = idf + "_" + tree_type;

		int region_id = data.reach_ids[rightOrder];
		boundaryTableRight.open(CTOutputLocation + "BoundaryTables/" + to_string(region_id) + "_" + idf + ".csv");
		boundaryTableRight << "SourceId" << "," << "Cost" << endl;

		float maxCost = MAXCOST;
		//std::cout << endl << data.reach_ids[rightOrder] << "->" << data.rightbfsOrder[rightOrder].size() << "->"; // TODO: UNCOMMENT
		int bfsroot = data.rightbfsRootNodes[rightOrder];
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}


			int nodeCls_new = nodeCls;
			/*if (data.allNodes[bfsroot]->isNa == 1) {
				nodeCls_new = -1;
			}*/

			/*if (data.allNodes[bfsroot]->p < 0.5) {
				nodeCls_new = 0;
			}*/

			/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls_new; // NC

			updatedNodes.push_back(data.allNodes[bfsroot]->originalId);

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}

	

					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[rightOrder] * 2;
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					// CHECK: set NA regions as -1
					int nodeCls_new_neigh = nodeCls;
					/*if (data.allNodes[neigh_id]->isNa == 1) {
						nodeCls_new_neigh = -1;
					}*/

					/*if (data.allNodes[neigh_id]->p < 0.5) {
						nodeCls_new_neigh = 0;
					}*/

					/*mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;*/
					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls_new_neigh; // NC

					updatedNodes.push_back(data.allNodes[neigh_id]->originalId);

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}


		vector<int> leafnodes;
		//step 1: get list of leave nodes
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int nodeId = data.rightbfsOrder[rightOrder][i];
			if (data.allNodes[nodeId]->childrenID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
		}
		//step 2: get chain length
		for (int i = 0; i < leafnodes.size(); i++) {
			int chainLength = 0;
			int nodeId = leafnodes[i];
			int leafnodeId = nodeId;
			int roughness = 1;
			float boundary_cost = -1.0;
			//float chainMaxCost = data.allNodes[nodeId]->fel; // NC: cost if elevation for NC
			float chainMaxCost = data.allNodes[nodeId]->cost; // NC: cost if elevation for NC
			pair<float, float> temp_maxCost_boundaryCost_pair;/* = make_pair(chainMaxCost, -1.0);*/
			pair<int, float> temp_cl_cost_pair /*= make_pair(0, -1.0)*/;
			pair<int, int> temp_chainLength_id/* = make_pair(-1, -1)*/;
			pair<float, int> temp_chainMaxCost_id/* = make_pair(chainMaxCost, -1)*/;
			pair<float, int> temp_cost_boundaryNode_pair /*= make_pair(-1.0,-1)*/;
			pair<int, vector<int>> temp_boundaryNode_leafNodes_pair;
			if (EB == 0) {
				temp_cl_cost_pair = make_pair(0, -1.0);
				temp_chainLength_id = make_pair(-1, -1);
			}
			else if (EB == 1) {
				temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);
				temp_chainMaxCost_id = make_pair(chainMaxCost, -1);
			}
			else if (EB == 2) {
				temp_cost_boundaryNode_pair = make_pair(-1.0, -1);
			}
			//pair<float, float> temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);

			while (data.allNodes[nodeId]->parentsID.size() != 0) {
				/*if (nodeClass[nodeId] == 1 && boundary_cost < 0.0) {*/
				if (nodeClass[nodeId] == 1 && boundary_cost < 0.0 && data.allNodes[nodeId]->isNa == 0) { // NC; do not consider NA nodes(nodes outside of RGB region)
					bool allchildObserved = true;
					bool observed = true;

					if (BOUNDARY_NODES_OBSERVED == 2) {  // uncomment after adding pits and Na identifiers
						for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
							int childid = data.allNodes[nodeId]->childrenID[idx];
							/*if (data.allNodes[childid]->isPits == 1 || data.allNodes[nodeId]->isNa == 1) {*/
							if (data.allNodes[childid]->isNa == 1) {
								allchildObserved = false;
								break;
							}
						}

						/*if (data.allNodes[nodeId]->isNa == 1 || data.allNodes[nodeId]->isPits ==1) {*/
						if (data.allNodes[nodeId]->isNa == 1) { // NC
							observed = false;
						}
					}
					/*if (data.allNodes[nodeId]->roughness == 1 || !observed || !allchildObserved || data.allNodes[nodeId]->isNa == 1) {
						break;
					}*/
					if (!observed || !allchildObserved || data.allNodes[nodeId]->isNa == 1) { // NC
						break;
					}
					boundary_cost = data.allNodes[nodeId]->cost;
					if (EB == 0) {
						temp_cl_cost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainLength_id.second = leafnodeId;
					}
					else if (EB == 1) {
						temp_maxCost_boundaryCost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainMaxCost_id.second = leafnodeId;
					}
					else if (EB == 2) {
						temp_cost_boundaryNode_pair.first = data.allNodes[nodeId]->cost;
						temp_cost_boundaryNode_pair.second = nodeId;

						if (data.boundaryNode_leafNodes_Map.find(nodeId) == data.boundaryNode_leafNodes_Map.end()) {
							vector<int> leafs;
							leafs.push_back(leafnodeId);
							data.boundaryNode_leafNodes_Map.insert(make_pair(nodeId, leafs));
						}
						else {
							data.boundaryNode_leafNodes_Map[nodeId].push_back(leafnodeId);
						}
					}
					data.leafNode_boundaryNodes.insert(make_pair(leafnodeId, nodeId));
					roughness = 0;
					if ((EB == 1) || (EB == 2)) {
						break;
					}
				}
				if (EB == 0) {
					chainLength++;
				}
				nodeId = data.allNodes[nodeId]->parentsID[0];
			}
			if (EB == 0) {
				temp_cl_cost_pair.first = chainLength;
				temp_chainLength_id.first = chainLength;
			}
			if (roughness == 0) {
				if (EB == 0) {
					data.chainLength_cost_Pairs.push_back(temp_cl_cost_pair);
					data.chainLength_nodeid_Pairs.push_back(temp_chainLength_id);
				}
				else if (EB == 1) {
					data.maxChainCost_cost_Pairs.push_back(temp_maxCost_boundaryCost_pair);
					data.maxCost_nodeid_Pairs.push_back(temp_chainMaxCost_id);
				}
				else if (EB == 2) {
					data.cost_boundaryNode_Pairs.insert(temp_cost_boundaryNode_pair);
				}
			}
		}
		if (EB == 0) {
			if (data.chainLength_cost_Pairs.size() != data.chainLength_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}
		else if (EB == 1) {
			if (data.maxChainCost_cost_Pairs.size() != data.maxCost_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}

		if (EB == 0) {
			// using chain length or max cost
			if (data.chainLength_cost_Pairs.size() != 0) {
				sort(data.chainLength_cost_Pairs.rbegin(), data.chainLength_cost_Pairs.rend());
				sort(data.chainLength_nodeid_Pairs.rbegin(), data.chainLength_nodeid_Pairs.rend());

				//top 20 percent
				//int top = data.chainLength_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.chainLength_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.chainLength_cost_Pairs[j].second;
					if (data.chainLength_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for chain length
					data.boundary_LeafNodes.push_back(data.chainLength_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostRight[rightOrder] = avg;
				}

				std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
			}
		}

		else if (EB == 1) {
			// using max cost branches
			if (data.maxChainCost_cost_Pairs.size() != 0) {

				sort(data.maxChainCost_cost_Pairs.rbegin(), data.maxChainCost_cost_Pairs.rend());
				sort(data.maxCost_nodeid_Pairs.rbegin(), data.maxCost_nodeid_Pairs.rend());


				//top 20 percent
				//int top = data.maxChainCost_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.maxChainCost_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.maxChainCost_cost_Pairs[j].second;
					if (data.maxCost_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for max cost
					data.boundary_LeafNodes.push_back(data.maxCost_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				/*if (avg > parameter.minCost) {
					data.inferedmaxCostRight[rightOrder] = avg;
				}*/

				data.inferedmaxCostRight[rightOrder] = avg; // NC

				if (EB == 1) {
					std::cout << data.maxChainCost_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
				}
				else {
					std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
				}

			}
		}

		else if (EB == 2) {
			//using boundary cost
			//
			// maintain vector of costs
			vector<float> boundaryCostList;

			if (data.cost_boundaryNode_Pairs.size() != 0) {
				if (DEBUG_OUTPUT == 1) {
					vector<float> infered_cost;
					set<pair<float, int>>::reverse_iterator it;
					for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
						pair<float, int> p = *it;
						infered_cost.push_back(p.first);
					}
				}

				// TODO: don't use cutoff for experiment
				//int top = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
				//int top = data.cost_boundaryNode_Pairs.size();

				int top = data.cost_boundaryNode_Pairs.size();
				int top_20 = 0;
				int top_80 = 0;
				if (parameter.useCutoff == 1) {
					top_20 = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
					top_80 = data.cost_boundaryNode_Pairs.size() * (1 - parameter.cutoff_percentage);
					top = top_80 - top_20;
				}

				float sum = 0.0;
				set<pair<float, int>>::reverse_iterator it;
				int counter = 0;
				for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
					// throw away top 20 and bottom 20 percent
					if (parameter.useCutoff == 1) {
						counter++;
						if (counter <= top_20) continue;
						if (counter > top_80) break;
					}

					pair<float, int> p = *it;
					int bNode = p.second;
					sum = sum + p.first; 

					boundaryCostList.push_back(p.first);

					// get children
					string children_id = "";
					for (int k = 0; k < data.allNodes[bNode]->childrenID.size(); k++) {
						int child_id = data.allNodes[bNode]->childrenID[k];
						children_id = children_id + "," + to_string(child_id);
					}

					// dump boundary cost to text file
					//boundary_costs_right << bNode << "," << p.first << "#" << data.allNodes[bNode]->parentsID[0] << "!" << children_id << "@";

					// dump boundary cost of each region to separate csv file
					boundaryTableRight << bNode << "," << p.first << endl;

					/*for (int n = 0; n < data.boundaryNode_leafNodes_Map[bNode].size(); n++) {
						int lNode = data.boundaryNode_leafNodes_Map[bNode][n];
						data.boundary_LeafNodes.push_back(lNode);
					}*/

					// don't use cutoff
					if (parameter.useCutoff != 1) {
						counter++;
						if (counter == top) {
							break;
						}
					}
				}
				float avg = -1.0;
				float q3 = -1.0;
				float stdev = -1.0;
				if (top != 0) {
					avg = sum / top;

					stdev = calc_std(boundaryCostList, sum, top);
					q3 = calc_q3(boundaryCostList);
					cout << "q3: " << q3 << endl;
					cout << "avg: " << avg << endl;

					//cout << "Reach ID: " << data.reach_ids[rightOrder] << " std: " << stdev << endl;

					//extra.standardDeviationRight[rightOrder] = stdev;

					//if (stdev > 2.0) {
					//	data.inferedmaxCostLeft[rightOrder] = -1;
					//	for (int mm = 0; mm < updatedNodes.size(); mm++) {
					//		mappredictions[updatedNodes[mm]] = -1; // reset updated predn to -1 // TODO: check
					//	}
					//}
				}


				//if (avg > parameter.minCost) {
				data.inferedmaxCostRight[rightOrder] = avg;
				//data.inferedmaxCostRight[rightOrder] = q3; // check with q3
				//} // NC
				std::cout << "->" << data.inferedmaxCostRight[rightOrder] << endl;

			}
		}



		data.cost_boundaryNode_Pairs.clear();
		data.boundaryNode_leafNodes_Map.clear();

		//boundary_costs_right << endl;
		boundaryTableRight.close();

	}
	//boundary_costs_right.close();
	cout << "Rightbank inference finished" << endl;
}


void cFlood::verify_deltaResult_left() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Leftbank: " << endl;


	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {

		vector<float> maxcost_accumulator;
		/*if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}*/

		// add reach ids of inferred regions
		data.leftInferredRegionsOld.push_back(data.reach_ids[leftOrder]);

		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[leftOrder] << "->" << data.leftbfsOrder[leftOrder].size() << "->"; // TODO: UNCOMMENT
		int bfsroot = data.leftbfsRootNodes[leftOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			/*if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				cout << "nodecls 1 for reachid: " << data.reach_ids[leftOrder] << endl;
				nodeCls = 1;
			}*/

			if (data.allNodes[bfsroot]->p > 0.5) {
				nodeCls = 1;
			}
			else {
				nodeCls = 0;
			}


			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					/*if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}*/

					if (data.allNodes[node]->p > 0.5) {
						cClass = 1;
					}
					else {
						cClass = 0;
					}

					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[leftOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}

			}
		}
	}
}

void cFlood::verify_deltaResult_right() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees
	std::cout << "Rightbank: " << endl;


	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {

		vector<float> maxcost_accumulator;
		/*if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}*/

		// add reach ids of inferred regions
		data.rightInferredRegions.insert(make_pair(data.reach_ids[rightOrder], true));

		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[rightOrder] << "->" << data.rightbfsOrder[rightOrder].size() << "->"; // TODO: UNCOMMENT
		int bfsroot = data.rightbfsRootNodes[rightOrder]; //// get the root
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			/*if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				cout << "nodecls 1 for reachid: " << data.reach_ids[leftOrder] << endl;
				nodeCls = 1;
			}*/

			if (data.allNodes[bfsroot]->p > 0.5) {
				nodeCls = 1;
			}
			else {
				nodeCls = 0;
			}

			//// TODO: added this to check result before inundation
			//if (nodeCls == 1) {
			//	mappredictions[data.allNodes[bfsroot]->originalId] = data.reach_ids[leftOrder];
			//}
			//else {
			//	mappredictions[data.allNodes[bfsroot]->originalId] = 0;
			//}

			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;

			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					/*if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}*/

					if (data.allNodes[node]->p > 0.5) {
						cClass = 1;
					}
					else {
						cClass = 0;
					}

					nodeClass[neigh_id] = cClass;

					/*if (nodeCls == 1) {
						mappredictions[data.allNodes[neigh_id]->originalId] = data.reach_ids[leftOrder];
					}
					else {
						mappredictions[data.allNodes[neigh_id]->originalId] = 0;
					}*/

					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;

					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}

			}
		}
	}
}

void cFlood::output() {
	auto start = std::chrono::system_clock::now();

	ofstream classout;
	string idf = "no_cutoff";
	if (parameter.useCutoff == 1) {
		idf = "cutoff_" + to_string((int)(parameter.cutoff_percentage * 100));
	}

	string tree_type = "fist_tree";
	if (parameter.useHMT == 1) {
		tree_type = "hmt_tree";
	}

	idf = idf + "_" + tree_type;

	classout.open(CTOutputLocation + "lambda_" + parameter.lambda_str.str() + "_" + CTPredictionTxt);
	//classout.open(CTOutputLocation + parameter.fname + ".txt");

	for (int i = 0; i < mappredictions.size(); i++) {
		classout << mappredictions[i] << endl;
	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
			prediction[row][col] = mappredictions[index];
			index++;

		}
	}

}



void cFlood::learning() {
	infer.marginal_ZnZpn.resize(parameter.allPixelSize * cNum * cNum); // All except bottom nodes
	infer.marginal_Zn.resize(parameter.allPixelSize * cNum); // Marginal Zn

	int iterateTimes = 0;
	bool iterator = true;

	double PiOld, EpsilonOld;
	double MuOld[cNum][Dim], SigmaOld[cNum][Dim][Dim];

	while (iterator) {
		this->UpdateTransProb();

		//copy current parameters to compare across iterations
		PiOld = parameter.Pi;
		EpsilonOld = parameter.Epsilon;
		for (int c = 0; c < cNum; c++) {
			for (int i = 0; i < Dim; i++) {
				MuOld[c][i] = parameter.Mu[c][i];
				for (size_t j = 0; j < Dim; j++) {
					SigmaOld[c][i][j] = parameter.Sigma[c][i][j];
				}
			}
		}

		cout << "Before MessagePropagation" << endl;
		this->MessagePropagation();

		cout << "Before UpdateMarginalProb" << endl;
		this->UpdateMarginalProb();

		cout << "Before UpdateParameters" << endl;
		this->UpdateParameters();

		cout << "Before UpdatePX_Z" << endl;
		this->UpdatePX_Z();

		//inference();

		//clock_t stop_s = clock();
		std::cout << endl << endl << "Iterate: " << iterateTimes /*<< "  Time: " << (stop_s - start_s) / float(CLOCKS_PER_SEC)*/;
		std::cout << endl << "Epsilon: " << eexp(parameter.Epsilon) << "  Pi: " << eexp(parameter.Pi);
		for (int c = 0; c < cNum; c++) {
			std::cout << endl << "Mu" << c;
			for (size_t i = 0; i < Dim; i++) {
				cout << " " << parameter.Mu[c][i] << " ";
			}
		}
		for (int c = 0; c < cNum; c++) {
			cout << endl << "Sigma" << c;
			for (size_t i = 0; i < Dim; i++) {
				for (size_t j = 0; j < Dim; j++) {
					cout << " " << parameter.Sigma[c][i][j] << " ";
				}
			}
		}

		//check stop criteria
		{
			bool MuConverge = true, SigmaConverge = true;
			double thresh = parameter.THRESHOLD;

			for (int c = 0; c < cNum; c++) {
				for (int i = 0; i < Dim; i++) {
					if (fabs((parameter.Mu[c][i] - MuOld[c][i]) / MuOld[c][i]) > thresh) {
						MuConverge = false;
						break;
					}

					for (int j = 0; j < Dim; j++) {
						if (fabs((parameter.Sigma[c][i][j] - SigmaOld[c][i][j]) / SigmaOld[c][i][j]) > thresh) {
							SigmaConverge = false;
							break;
						}
					}
				}
			}

			double epsilonRatio = fabs((eexp(parameter.Epsilon) - eexp(EpsilonOld)) / eexp(EpsilonOld));
			double PiRatio = fabs((eexp(parameter.Pi) - eexp(PiOld)) / eexp(PiOld));
			//double MRatio = fabs((eexp(parameter.M) - eexp(MOld)) / eexp(MOld));

			if (epsilonRatio < thresh && PiRatio < thresh && MuConverge && SigmaConverge) {
				iterator = false;
			}

			iterateTimes++;
			//cout << "Iteration " << iterateTimes << " Done.." << endl;
			if (iterateTimes >= parameter.maxIteratTimes) {
				iterator = false;
			}
		}

	} // end while

}

int main(int argc, char* argv[]) {
	cFlood flood;
	flood.input(argc, argv);
	return 0;
}


//Testing Functions

void cFlood::getOriginIdBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg * 2;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getIds() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = nodid;
			}
		}
	}

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = nodid;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_Ids.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOrgIds() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = data.allNodes[nodid]->originalId;
			}
		}
	}

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = data.allNodes[nodid]->originalId;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_OrgIds.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOriginIdLeftBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdLeftBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOriginIdRightBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdRightBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int regionSize = data.leftbfsOrder[leftOrder].size();
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int regionSize = data.rightbfsOrder[rightOrder].size();
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_regionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getLeftRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int regionSize = data.leftbfsOrder[leftOrder].size();
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}

	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_LeftregionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}

void cFlood::getRightRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int regionSize = data.rightbfsOrder[rightOrder].size();
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_RightregionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::sanityChecker() {
	//for left trees
	//std::cout << "Leftbank: " << endl;
	std::cout << "Starting Sanity Checks..." << endl;
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}
		// first find list of leaves.
		vector<int> leavesList;
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int curr = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[curr]->childrenID.size() == 0) {
				leavesList.push_back(curr);
			}
		}
		for (int i = 0; i < leavesList.size(); i++) {
			int curr = leavesList[i];
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] == 1) {
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
			}
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] != 1) {
					std::cout << i << " !!!Error Elevation Assumption Failed!!!!" << endl;
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
				//std::cout << "curr= " << curr<<endl;
			}
		}
	}

	//for right trees
	//std::cout << "Rightbank: " << endl;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		// first find list of leaves.
		vector<int> leavesList;
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int curr = data.rightbfsOrder[rightOrder][i];
			//std::cout << "curr = " << curr << endl;
			if (data.allNodes[curr]->childrenID.size() == 0) {
				leavesList.push_back(curr);
			}
		}
		for (int i = 0; i < leavesList.size(); i++) {
			int curr = leavesList[i];
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] == 1) {
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
			}
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] != 1) {
					std::cout << i << " !!!Error Elevation Assumption Failed!!!!" << endl;
					break;
				}
				else {
					curr = data.allNodes[curr]->parentsID[0];
				}
			}
		}
	}
	std::cout << "Sanity Test Complete..." << endl;
}

void cFlood::getOriginIdBanks_effectiveBranches() {
	cout << "inside getOriginIdBanks_effectiveBranches function" << endl;
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}

	//for right trees
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg * 2;
			}
		}
	}

	for (int i = 0; i < data.boundary_LeafNodes.size(); i++) {
		int nodeId = data.boundary_LeafNodes[i];
		int leafNode = nodeId;
		tempMap[data.allNodes[nodeId]->originalId] = 2;
		nodeId = data.allNodes[nodeId]->parentsID[0];
		//std::cout << "i= " << i << " nodeId = " << nodeId << endl;
		while (data.allNodes[nodeId]->parentsID.size() != 0) {
			tempMap[data.allNodes[nodeId]->originalId] = 1;
			nodeId = data.allNodes[nodeId]->parentsID[0];
		}
		if (data.leafNode_boundaryNodes.find(leafNode) != data.leafNode_boundaryNodes.end()) {
			int bNode = data.leafNode_boundaryNodes[leafNode];
			tempMap[data.allNodes[bNode]->originalId] = 3;
		}
	}
	cout << "getOriginIdBanks_effectiveBranches function ended" << endl;

	//ofstream classout;
	//classout.open(CTOutputLocation + parameter.fname + "_Viz2.txt");
	//for (int i = 0; i < tempMap.size(); i++) {
	//	classout << tempMap[i] << endl;
	//}
	//classout.close();
}



// Input:	waterlikelihood: the probabilty for water node from Unet.(Not loglikelihood) The order should follow the tree sturcture.
//			drylikelihood: the probabilty for water node from Unet.The order should follow the tree sturcture.
//			treeInput: the Id of the node in tree structure. The are only nodes in main branch. The order of the nodes is from lower to higher(i.e. the first Id should be the id of the lowest node in one region, the river node).
//			Nodes: The information for parents and children.
//			treeLength: the length of the tree structure for traversing the tree structure. 
// Output:	The vector of the loglikelihood for every frontier nodes in one region.
// Global parameter: rho and Pi. In code, the names are parameter.Pi and parameter.rho. Please modify the names of these parameters
void cFlood::getLoglikelihood() {
	parameter.Pi = parameter.Pi_orig;
	parameter.rho = 0.999;
	//get Gain
	double initialLog = 0;




	// Right Bank
	int rightOrder_selected = 0;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		int numInt = 0;
		double curGain = 0;
		double curWaterProb, curDryProb, curMaxGain = 0;
		vector<double> loglikelihood;
		vector<int> mainBranchNodeIds;
		vector<Node*> sortedbfsOrder;
		for (int j = 0; j < data.rightbfsOrder[rightOrder].size(); j++)
		{
			sortedbfsOrder.push_back(data.allNodes.at(data.rightbfsOrder[rightOrder][j]));
		}
		std::sort(sortedbfsOrder.begin(), sortedbfsOrder.end(), comp);
		for (int i = 0; i < sortedbfsOrder.size(); i++)
		{
			bool parentsAllVisited = true;
			double curWaterProb, curDryProb, curMaxGain = 0;

			Node* currNode;
			currNode = sortedbfsOrder[i];
			currNode->visited = true;
			parentsAllVisited = true;
			if (currNode->isNa == 1) {
				curDryProb = eln_ll(0.5);
				curWaterProb = eln_ll(0.5);
			}
			else {
				curDryProb = eln_ll(1 - currNode->p);
				curWaterProb = eln_ll(currNode->p);
			}
			//cout << curGain << ',' << curDryProb << ',' << eln_ll(parameter.Pi) << endl;
			if (currNode->parentsID.size() == 0) {
				if (currNode->childrenID.size() != 0) {
					if (data.allNodes.at(currNode->childrenID[0])->parentsID.size() == 1)
						currNode->curGain = curWaterProb - curDryProb - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho); // TODO: check this +/-
					else {
						for (size_t a = 0; a < data.allNodes.at(currNode->childrenID[0])->parentsID.size(); a++) {
							/*if (data.allNodes[data.allNodes.at(currNode->childrenID[0])->parentsID[a]]->visited == true) {
								if (data.allNodes.at(currNode->childrenID[0])->parentsID[a] != currNode->nodeIndex)
								{
									parentsAllVisited = false;
									break;
								}
							}*/
							if (data.allNodes[data.allNodes.at(currNode->childrenID[0])->parentsID[a]]->visited == false) {
								parentsAllVisited = false;
								break;
							}
						}
						if (parentsAllVisited == true) {
							currNode->curGain = curWaterProb - curDryProb - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho); // TODO: check this +/-
						}
						//currNode->curGain = curWaterProb - curDryProb - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho); // TODO: check this +/-
						else {
							currNode->curGain = curWaterProb - curDryProb - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi);
						}
					}

				}
				else {//false node which is seperated from other nodes, chain length = 1, may due to mosaic or data error
					currNode->curGain = 0;
				}
			}
			else {
				if (currNode->childrenID.size() != 0) {
					if (data.allNodes.at(currNode->childrenID[0])->parentsID.size() == 1)
					{
						currNode->curGain = curWaterProb - curDryProb + eln_ll(parameter.rho) - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho);
						//currNode->curGain = curWaterProb - curDryProb + parameter.rho - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho);
					} // TODO: check this +/-
					else {
						for (size_t a = 0; a < data.allNodes.at(currNode->childrenID[0])->parentsID.size(); a++) {
							/*if (data.allNodes[data.allNodes.at(currNode->childrenID[0])->parentsID[a]]->visited == true) {
								if (data.allNodes.at(currNode->childrenID[0])->parentsID[a] != currNode->nodeIndex)
								{
									parentsAllVisited = false;
									break;
								}
							}*/
							if (data.allNodes[data.allNodes.at(currNode->childrenID[0])->parentsID[a]]->visited == false) {
								parentsAllVisited = false;
								break;
							}
						}
						if (parentsAllVisited == true) {
							currNode->curGain = curWaterProb - curDryProb + eln_ll(parameter.rho) - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho);
							//currNode->curGain = curWaterProb - curDryProb + parameter.rho - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi);

						}
						else {
							currNode->curGain = curWaterProb - curDryProb + eln_ll(parameter.rho) - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi);
						}
						//currNode->curGain = curWaterProb - curDryProb + parameter.rho - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi) + eln_ll(1 - parameter.rho);
					}
				}
				else
					currNode->curGain = curWaterProb - curDryProb + eln_ll(parameter.rho) - eln_ll(1 - parameter.rho) - eln_ll(parameter.Pi) + eln_ll(1 - parameter.Pi);
			}
		}
		for (int i = 0; i < sortedbfsOrder.size(); i++)
		{
			sortedbfsOrder[i]->visited = false;
		}
		//int parentId = curNode;
		initialLog = 0;
		for (int j = 0; j < data.rightbfsOrder[rightOrder].size(); j++) {
			if (data.allNodes[data.rightbfsOrder[rightOrder][j]]->isNa == 1) {
				curDryProb = eln_ll(0.5);
				curWaterProb = eln_ll(0.5);
			}
			else {
				curDryProb = eln_ll(1 - data.allNodes[data.rightbfsOrder[rightOrder][j]]->p);
				curWaterProb = eln_ll(data.allNodes[data.rightbfsOrder[rightOrder][j]]->p);
			}

			curGain = curDryProb - eln_ll(1 - parameter.Pi);

			initialLog += curGain;
		}


		int curNode, curNodePosition;
		curNode = data.rightNodesInOrder[rightOrder];
		curNodePosition = 0;
		int parentId = curNode;
		queue<int> numQue;
		curMaxGain = initialLog;
		int flag = 0;
		int numNode = 0;
		while (data.allNodes[curNode]->childrenID.size() != 0) {
			if (curNodePosition < sortedbfsOrder.size())
				while (sortedbfsOrder[curNodePosition]->cost < data.allNodes[curNode]->cost || sortedbfsOrder[curNodePosition]->cost == data.allNodes[curNode]->cost)
					//while (sortedbfsOrder[curNodePosition]->cost < data.allNodes[curNode]->cost)
				{
					sortedbfsOrder[curNodePosition]->visited = true;
					numQue.push(curNodePosition);
					curNodePosition++;
					/*numNode++;
					if (numNode > 0)
						break;*/
					if (curNodePosition < sortedbfsOrder.size())
						continue;
					else
						break;
				}
			/*else
			{
				cout << curNode << endl;
			}*/

			//if (flag == 0)
			//{
			//	/*for (int a = 0; a < sortedbfsOrder.size(); a++)
			//	{
			//		for (int b = 0; b < sortedbfsOrder[a]->parentsID.size(); b++)
			//			if (sortedbfsOrder[a]->parentsID.size()!=0)
			//				if (sortedbfsOrder[a]->cost < data.allNodes[sortedbfsOrder[a]->parentsID[b]]->cost)
			//				{
			//					cout << a << endl;
			//				}
			//	}*/
			//	//cout << data.allNodes[sortedbfsOrder[curNodePosition]->parentsID[0]]->cost - data.allNodes[data.leftNodesInOrder[leftOrder]]->cost << endl;
			//	cout << sortedbfsOrder[curNodePosition - 1]->parentsID.size() << endl;
			//	cout << sortedbfsOrder[curNodePosition - 1]->childrenID.size() << endl;
			//	flag = 1;
			//}

			while (!numQue.empty())
			{
				int curNodeNum = numQue.front();
				numQue.pop();

				curMaxGain = curMaxGain + sortedbfsOrder[curNodeNum]->curGain;
				curGain = curDryProb - eln_ll(1 - parameter.Pi);

			}
			loglikelihood.push_back(curMaxGain);
			mainBranchNodeIds.push_back(curNode);
			parentId = curNode;
			curNode = data.allNodes[curNode]->childrenID[0];

		}
		data.loglikelihood_rightRegions.push_back(loglikelihood);
		data.mainBranchNodeIds_rightRegions.push_back(mainBranchNodeIds);
		data.rightOrderSelectedToRightOrder[rightOrder_selected] = rightOrder;
		rightOrder_selected++;
	}
}

void cFlood::viterbi(double lambda, int ranges) {
	vector<double> llCurrRegion;
	vector<double> llNextRegion;
	vector<double> costMapCurrRegion;
	vector<double> costMapNextRegion;
	// 	vector<vector<double>> delta;
	vector<vector<pair<double, int>>> delta;
	vector<vector<pair<double, int>>> delta1;
	vector<vector<int>> path;
	vector<map<int, int>> path_map;
	// 	vector<double>* tempLoss;
		// vector<double> tempLoss;
	vector<pair<double, int>> tempLoss;
	double loss, max_val;
	double tempBest;
	int individualPath, sta;
	// vector<size_t> result_ids_left;
	// vector<double> result_costs;

	// TODO: check this parameter
	// parameter.lambda = 0.5;
	size_t region_size = data.loglikelihood_leftRegions.size();
	//size_t region_size = 10;
	size_t rangeSize;
	size_t rangeSizeN;
	ofstream filename,filename1;
	double prev_lb, prev_ub = 0;
	prev_lb = 0;
	vector<int> regularizedSizeLeft, regularizedSizeRight;
	size_t bestRange[] = { 0,0 };
	size_t rangeNum = 2;
	vector<int> lb, ub;
	auto start = std::chrono::system_clock::now();
	// for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++)
	// {
	// 	if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
	// 		continue;
	// 	}
	// 	regularizedSizeLeft.push_back(data.leftbfsOrder[leftOrder].size());
	// }
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++)
	{
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		regularizedSizeRight.push_back(data.rightbfsOrder[rightOrder].size());
	}


	//Right Region
	region_size = data.loglikelihood_rightRegions.size();
		for (int rightOrder = 0; rightOrder < region_size-1; rightOrder++) {
			vector<pair<double, int>> temp;
			vector<int> tempPath;
			map<int, int> tempPathMap;
			llCurrRegion = data.loglikelihood_rightRegions[rightOrder];
			llNextRegion = data.loglikelihood_rightRegions[rightOrder + 1];
			rangeSize = llCurrRegion.size()/ ranges; //Set the size of range
			if (llCurrRegion.size()/ ranges < 1)
				rangeSize = 1;
			rangeSizeN =  llNextRegion.size() / ranges;
			if (llNextRegion.size()/ ranges < 1)
				rangeSizeN = 1;
			if (rightOrder != 0)
				tempLoss = delta.back();
			for (int j = 0; j < llNextRegion.size(); j++) {

				size_t n = j * rangeSizeN;

				if (n > llNextRegion.size() || n == llNextRegion.size())
				{
					//cout << m << endl;
					n = llNextRegion.size() - 1;
					//break;
				}
				double nextCost = data.allNodes[data.mainBranchNodeIds_rightRegions[rightOrder + 1][n]]->cost;
				tempBest = MAXGAIN;
				//cout << j << endl;




				for (int k = 0; k < llCurrRegion.size(); k++) {
					size_t m = k * rangeSize;

					if (m > llCurrRegion.size() || m == llCurrRegion.size())
					{
						//cout << m << endl;
						m = llCurrRegion.size() - 1;
						//break;
					}
					double currCost = data.allNodes[data.mainBranchNodeIds_rightRegions[rightOrder][m]]->cost;
					if (rightOrder == 0) {
						loss = -llCurrRegion[m] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost);
					}
					else if (rightOrder == region_size - 2)
					{

						loss = tempLoss[k].first - llCurrRegion[m] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost) - llNextRegion[n] / llNextRegion.size();
					}
					else {
						// 	tempLoss = &delta.back();
						loss = tempLoss[k].first - llCurrRegion[m] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost);
					}
					if (loss < tempBest) {
						tempBest = loss;
						individualPath = k; // save the best route and remove all other routes
					}
					if (m == llCurrRegion.size() - 1)
					{
						break;
					}
				}

				temp.push_back(make_pair(tempBest, j));
				tempPath.push_back(individualPath); // save the best route in one region
				tempPathMap[j] = individualPath;
				if (n == llNextRegion.size() - 1)
				{
					break;
				}
			}
// 			cout << tempPath.size() << endl;

// 			cout << "Before delta" << endl;
			delta.push_back(temp);
			path.push_back(tempPath);
			path_map.push_back(tempPathMap);
		}
// 		cout << "delta back size: " << delta.back().size() << endl;
// 		cout << "path size: " << path.size() << endl;

		// find one final optimal path
		max_val = delta.back().at(0).first;
		sta = 0;

// 		cout << "max_val: " << max_val << endl;

		// find frontier node (Start point) on last region
		for (int i = 0; i < delta.back().size(); i++) {
			if (delta.back().at(i).first < max_val) {
				max_val = delta.back().at(i).first;
				sta = i;
			}
		}
// 		cout << "max val: " << max_val << endl;
// 		cout << "sta: " << sta << endl;
		// 	data.result_ids_right.push_back(sta); // TODO: should not be 0
		llCurrRegion = data.loglikelihood_rightRegions[region_size - 1];
		rangeSize =  llCurrRegion.size()/ ranges;
		if (llCurrRegion.size() / ranges < 1)
		{
			rangeSize = 1;
		}
		ub.push_back(sta * rangeSize + rangeSize);
		if (ub.back() > llCurrRegion.size())
			ub.back() = llCurrRegion.size();

		lb.push_back(sta * rangeSize - rangeSize);
		//if (sta * rangeSize - rangeSize > llCurrRegion.size())
		//	lb.back() = llCurrRegion.size();
		if (lb.back() < 0)
			lb.back() = 0;
		//if (sta * rangeSize - rangeSize > sta * rangeSize + rangeSize)
		//	lb.back() = ub.back();
// 		cout << "sta corr: " << delta.back().at(sta).second << endl;
// 		cout << rangeSize << endl;
// 		cout << ub.back() << endl;
		//  	for (int rightOrder=0; rightOrder < data.rightNodesInCorrectOrder.size()-1; rightOrder++){
		for (int rightOrder = 0; rightOrder < region_size - 1; rightOrder++) {
			// 		sta = path[data.rightNodesInCorrectOrder.size() - rightOrder - 1][sta];
			llCurrRegion = data.loglikelihood_rightRegions[region_size - 1 - rightOrder - 1];
			if (rightOrder == 0) {
				sta = path[region_size - 1 - rightOrder - 1][sta];
			}
			else {
				sta = path_map[region_size - 1 - rightOrder - 1][sta];
			}
			// sta = path[2 - rightOrder - 1][sta];
// 			cout << "rightOrder: " << rightOrder << " sta: " << sta << " path size: " << path[region_size - 1 - rightOrder - 1].size() << endl;
			rangeSize = llCurrRegion.size() / ranges;

			if (llCurrRegion.size() / ranges < 1)
			{
				rangeSize = 1;
			}
			ub.push_back(sta * rangeSize + rangeSize);
// 			cout << rangeSize << endl;
// 			cout << ub.back() << endl;
			if (ub.back() > llCurrRegion.size())
				ub.back() = llCurrRegion.size();
			lb.push_back(sta * rangeSize - rangeSize);
			if (lb.back() < 0)
				lb.back() = 0;
		}
		reverse(lb.begin(), lb.end());
		reverse(ub.begin(), ub.end());
// 		for (int rightOrder = 0; rightOrder < lb.size(); rightOrder++)
// 			cout << lb[rightOrder] << "," << ub[rightOrder] << endl;
		path.clear();
		delta.clear();
		path_map.clear();
// 		cout << data.loglikelihood_rightRegions.size() << endl;
		for (int rightOrder = 0; rightOrder < region_size - 1; rightOrder++) {
			vector<pair<double, int>> temp;
			vector<int> tempPath;
			map<int, int> tempPathMap;
			llCurrRegion = data.loglikelihood_rightRegions[rightOrder];
			llNextRegion = data.loglikelihood_rightRegions[rightOrder + 1];
			if (rightOrder != 0)
				tempLoss = delta.back();
			for (int j = lb[rightOrder + 1]; j < ub[rightOrder + 1]; j++) {
				double nextCost = data.allNodes[data.mainBranchNodeIds_rightRegions[rightOrder + 1][j]]->cost;

				tempBest = MAXGAIN;
				//cout << j << endl;




				for (int k = lb[rightOrder]; k < ub[rightOrder]; k++) {
					double currCost = data.allNodes[data.mainBranchNodeIds_rightRegions[rightOrder][k]]->cost;



					if (rightOrder == 0) {
						loss = -llCurrRegion[k] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost);
					}
					else if (rightOrder == region_size - 2)
					{
						loss = tempLoss[k].first - llCurrRegion[k] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost) - llNextRegion[j] / llNextRegion.size();
					}
					else {
						// 	tempLoss = &delta.back();
						loss = tempLoss[k].first - llCurrRegion[k] / regularizedSizeRight[rightOrder] + lambda * abs(currCost - nextCost);
					}
					if (loss < tempBest) {
						tempBest = loss;
						individualPath = k; // save the best route and remove all other routes
					}
				}

				temp.push_back(make_pair(tempBest, j));
				tempPath.push_back(individualPath); // save the best route in one region
				tempPathMap[j] = individualPath;
			}
			//cout << tempPath.size() << endl;

// 			cout << "Before delta" << endl;
			delta.push_back(temp);
			path.push_back(tempPath);
			path_map.push_back(tempPathMap);
		}
// 		cout << "delta back size: " << delta.back().size() << endl;
// 		cout << "path size: " << path.size() << endl;

		// find one final optimal path
		max_val = delta.back().at(0).first;
		sta = 0;

// 		cout << "max_val: " << max_val << endl;

		// find frontier node (Start point) on last region
		for (int i = 0; i < delta.back().size(); i++) {
			if (delta.back().at(i).first < max_val) {
				max_val = delta.back().at(i).first;
				sta = i;
			}
		}
// 		cout << "max val: " << max_val << endl;
// 		cout << "sta: " << sta << endl;
		// 	data.result_ids_right.push_back(sta); // TODO: should not be 0
		data.result_ids_right.push_back(delta.back().at(sta).second);
// 		cout << "sta corr: " << delta.back().at(sta).second << endl;

		//  	for (int rightOrder=0; rightOrder < data.rightNodesInCorrectOrder.size()-1; rightOrder++){
		for (int rightOrder = 0; rightOrder < region_size - 1; rightOrder++) {
			// 		sta = path[data.rightNodesInCorrectOrder.size() - rightOrder - 1][sta];

			if (rightOrder == 0) {
				sta = path[region_size - 1 - rightOrder - 1][sta];
			}
			else {
				sta = path_map[region_size - 1 - rightOrder - 1][sta];
			}
			// sta = path[2 - rightOrder - 1][sta];
// 			cout << "rightOrder: " << rightOrder << " sta: " << sta << " path size: " << path[region_size - 1 - rightOrder - 1].size() << endl;
			data.result_ids_right.push_back(sta);
		}
		reverse(data.result_ids_right.begin(), data.result_ids_right.end());

		// display the result
		//  	for (int rOrder=0; rOrder < data.rightNodesInCorrectOrder.size(); rOrder++){
		for (int rOrder = 0; rOrder < region_size; rOrder++) {

			int region_id = data.rightNodesInCorrectOrder[rOrder];
			//cout << "Region Id: " << region_id << endl;
			//cout << "Number of Region Nodes:" << data.rightbfsOrder[data.rightOrderSelectedToRightOrder[rOrder]].size() << endl;
			//cout << "Max Elevation in this region:" << data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder].back()]->cost << endl;
			int frontier_node_idx = data.result_ids_right[rOrder];
			//cout << "Max Sum loglikelihood:" << data.loglikelihood_rightRegions[rOrder][numBranch] << " Max Sum Elevation:" << data.inferredFloodFrontier_regionId2Cost[data.rightNodesInCorrectOrder[rOrder]] << endl;
			double regularized_cost = data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder][frontier_node_idx]]->cost;
			//cout << " Viterbi Loglikelihood:" << data.loglikelihood_rightRegions[rOrder][frontier_node_idx] << " Viterbi Elevation: " << regularized_cost << endl;
			//int nextFrontierNode;
			//double maxLog = *max_element(data.loglikelihood_rightRegions[rOrder].begin(), data.loglikelihood_rightRegions[rOrder].end());
			//cout << "Maximum Loglikelihood:" << maxLog;
			//for (int i = 0; i < data.loglikelihood_rightRegions[rOrder].size(); i++)
			//{
			//	if (data.loglikelihood_rightRegions[rOrder][i] == maxLog)
			//	{
			//		cout << "Maximum Loglikelihood Elevation:" << data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder][i]]->cost << endl;
			//		break;
			//	}
			//}
			//if (rOrder == region_size - 1)
			//	nextFrontierNode = data.result_ids_right[rOrder];
			//else
			//	nextFrontierNode = data.result_ids_right[rOrder + 1];
			//if (rOrder == region_size - 1)
			//	cout << "RightMost" << endl;
			//else
			//{
			//	cout << "Elevation difference:" << abs(regularized_cost - data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder + 1][nextFrontierNode]]->cost) << "Normalization Value + Regularization term + :" << data.loglikelihood_rightRegions[rOrder][frontier_node_idx] / regularizedSizeRight[rOrder] << "+" << lambda * abs(regularized_cost - data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder + 1][nextFrontierNode]]->cost) << endl;
			//	cout << " Regularization term:" << lambda * abs(regularized_cost - data.allNodes[data.mainBranchNodeIds_rightRegions[rOrder + 1][nextFrontierNode]]->cost) << " Normalization Value:" << data.loglikelihood_rightRegions[rOrder][frontier_node_idx] / regularizedSizeRight[rOrder] << endl;
			//}
			data.regularizedMaxCostRight[data.rightOrderSelectedToRightOrder[rOrder]] = regularized_cost;
// 			cout << endl;
			// 		data.regularizedMaxCostright.push_back(regularized_cost);
		}
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double>elapsed_seconds = end - start;
		timeSave << elapsed_seconds.count() << endl;
		cout << "viterbi" << elapsed_seconds.count() << "Seconds" << endl;
// 		cout << ranges << endl;
// 		cout << endl;
		path.clear();
		path_map.clear();
		delta.clear();
		data.result_ids_right.clear();
		lb.clear();
		ub.clear();

}



void cFlood::saveExtraInfo() {
	ofstream idtable;
	idtable.open(CTOutputLocation + parameter.reachId + "_id_table.csv");

	idtable << "new_id" << "," << "original_id" << "," << "is_large_region" << "," << "is_river_id" << "," << "is_reach_id" << endl;

	for (int row = 0; row < parameter.ROW; row++) {
		for (int col = 0; col < parameter.COLUMN; col++) {
			int node_id_orig = row * parameter.COLUMN + col;
			int node_id_new = nodeIndexMap[node_id_orig];
			int is_large = 0;

			if (std::find(data.leftNodesInCorrectOrder.begin(), data.leftNodesInCorrectOrder.end(), node_id_new) != data.leftNodesInCorrectOrder.end() ||
				std::find(data.rightNodesInCorrectOrder.begin(), data.rightNodesInCorrectOrder.end(), node_id_new) != data.rightNodesInCorrectOrder.end()
				) {
				is_large = 1;
			}

			idtable << node_id_new << "," << node_id_orig << "," << is_large << "," << data.river_ids_map[node_id_orig] << "," << data.reach_ids_orig_map[node_id_orig] << endl;
		}
	}
	idtable.close();
}




void cFlood::clear_all() {
	cout << "Clearing all the vectors!";


}
