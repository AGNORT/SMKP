#ifndef LABELSETTING_H
#define LABELSETTING_H

#ifdef __cplusplus
extern "C" {
#endif
#include "cli_parse.h"
#include "instance.h"
#include "data_structs.h"
#include "item.h"
#include "ts_heur.h"
#include <string.h>
#ifdef __cplusplus
}
#endif

#define EX 1e-8
#define MAXMUMSOLTIME 3600

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <stack>
#include <bitset>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <thread>
#include <functional>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>
#include "gurobi_c++.h"

/*parameters*/
extern int g_tolItemNum;							//the total number of items
extern int g_firstDim;								//the first dimension of the bucket
extern int secondDim;								//the second dimension of the bucket
extern double g_bestLB;								//the best lower bound		
extern int g_secondDimDiv;							//the divisor for second dimension of the bucket
extern double g_smallest_bDiva;						//the smallest ratio of a/b for the present instance	
extern bool g_solveKnapsack;						//if solve the knapsack problem

extern Instance g_instance;							//g_instance data				

/*results*/
extern int g_tolGeneratedLabel;

/*bucket index*/
extern std::vector<int> g_item_to_firstIdx;			//map item to item bucket index
extern std::vector<int> g_weight_to_secondIdx;			//map weight to capacity bucket index
extern int g_domWidth;									//the width of buckets dominance check
extern bool g_pInteger;								//remark if parameter p is integer 
extern bool g_aInteger;								//remark if parameter a is integer 


typedef struct {
	size_t index;		//index of items in the orginstances
	size_t pos;			//index of deleted items in the g_instances
	double sum_p_ab;
	double profit;
} ItemIndex;

//The label class
class MyLabel {
public:
	MyLabel() {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = 0;
		parentLab = nullptr;
		//bitSet.resize(g_tolWrd);
	}
	/*MyLabel(int sr3Size) {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = 0;
		parentLab = nullptr;
		sr3Vis.resize(sr3Size);
	}*/
	MyLabel(double& binNumDual) {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = binNumDual;
		parentLab = nullptr;
		//bitSet.resize(g_tolWrd);
	}
	MyLabel(const MyLabel* preLab) {
		lastItem = preLab->lastItem;
		tolWeight = preLab->tolWeight;
		sum_a = preLab->sum_a;
		sum_b = preLab->sum_b;
		tolProfit = preLab->tolProfit;
		parentLab = preLab->parentLab;
		itemSet = preLab->itemSet;
		//sr3Vis = preLab->sr3Vis;
		//bitSet = preLab->bitSet;
	}

	//compare two vectors with elements only take 0 or 1
	bool areEqual(const std::vector<int>& a, const std::vector<int>& b) {
		std::bitset<1024> ba, bb;
		for (size_t i = 0; i < a.size(); ++i) {
			ba[i] = a[i];
			bb[i] = b[i];
		}
		return ba == bb;
	}

	bool operator==(const MyLabel& other) {
		return areEqual(this->itemSet, other.itemSet);
	}
public:
	int lastItem;           //the last inserted item in the bag
	double sum_a;			//the summation of a part
	double sum_b;			//the summation of b part
	double tolWeight;       //the total weight of the items in the bag
	double tolProfit;			//the total profit of the items in the bag
	MyLabel* parentLab;		//the parent label of the current label
	std::vector<int> itemSet;	//the inserted items in the bag, Size = g_tolItemNum, value = 1 if a item is inserted
	//std::vector<int> sr3Vis;  //the visiting time to each sr3
	//vector<int> bitSet;     //the bit projection of the inserted items in the bag
};
using Bin = MyLabel;		// rename MyLabel as Bin

//the best kanpsack solution
class KnapsackSol {
public:

public:
	MyLabel* bestLab;			//the best label
	std::string bestItemSet;	    //the best item set
};


/*function decleration*/
//convert vector to string
std::string JoinVector(const std::vector<int>& vec);

// instance data������,�Ӵ�С����
void sort_instance_by_p_ab(Instance& preInstances, ItemIndex*& indices);

//design the labelsetting algorithm to solve the submodular knapsack problem
int LabelSettingSolveKnapsack(Args& args);

//label extention
void LabelExtention(Instance& preInstances, MyLabel* parentLab, MyLabel* preLab, int item);

//completion bound
int CompletionBound(
	MyLabel* lab,
	int prefirstIdx,
	DblMatrix* ub_matr
);

//The whole dominance logic
bool DominanceLogic(
	int preSecondIdx,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>*& newExtended,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>*& oldExtended,
	MyLabel* preLab,
	bool newLabelFalg,
	bool heuDom);

void JgeDominance(
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>& newExtended,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>& oldExtended,
	MyLabel* lab,
	std::vector<MyLabel*>& dominatedOldLabs,
	bool newLabelFlag,
	bool heuDom,
	bool dualCalFlag);

int FKP(const Instance* RESTRICT const inst,
	DblMatrix* RESTRICT const output,
	const std::vector<double>& weightRec,
	bool dualFlag);


#endif
