/****************The second version****************/
#include "labelSetting.h"
#include "compactModel.h"

using namespace std;

Instance g_instance = { 0 };					//g_instance data

/*parameters*/
int g_tolItemNum = 0;							//the total number of items
int g_firstDim = 0;								//the first dimension of the bucket
int secondDim = 0;								//the second dimension of the bucket
double g_bestLB = 0;								//the best lower bound		
int g_secondDimDiv = 1;							//the divisor for second dimension of the bucket
double g_smallest_bDiva = 1e10;					//the smallest ratio of a/b for the present instance	
bool g_solveKnapsack = false;					//if solve the knapsack problem

vector<double> g_profitRec;						//the profit record from each item

/*bucket index*/
std::vector<int> g_item_to_firstIdx;			//map item to item bucket index
std::vector<int> g_weight_to_secondIdx;			//map weight to capacity bucket index
int g_domWidth = 6;								//the width of buckets dominance check
bool g_pInteger = true;						//remark if parameter p is integer 
bool g_aInteger = true;						//remark if parameter a is integer 

/*results*/
int g_tolGeneratedLabel = 0;					//the total number of generated labels
int g_nonDominatedLabel = 0;					//the total number of non-dominated labels
int g_CBFathomLabel = 0;						//the number of labels fathomed by CB
int g_dominatedLabel = 0;						//the number of labels fathomed by CB

/*sort the instances according to the sum of a and b*/
// ascending sort
int compare_item_index_ab(const void* a, const void* b) {
	double diff = ((ItemIndex*)b)->sum_p_ab - ((ItemIndex*)a)->sum_p_ab;
	return (diff > 0) - (diff < 0);
}
int compare_item_index_prt(const void* a, const void* b) {
	double diff = ((ItemIndex*)b)->profit - ((ItemIndex*)a)->profit;
	return (diff > 0) - (diff < 0);
}
// instance data
void sort_instance_by_p_ab(Instance& preInstances, ItemIndex*& indices) {
	size_t n = preInstances.n_items;
	indices = (ItemIndex*)malloc(n * sizeof(ItemIndex));
	
	for (size_t i = 0; i < n; ++i) {
		indices[i].index = i;
		indices[i].sum_p_ab = preInstances.p_ptr[i] / (preInstances.a_ptr[i] + sqrt(preInstances.b_ptr[i]));
		//indicesRec[i].sum_p_ab = 1 / (preInstances.a_ptr[i] + sqrt(preInstances.b_ptr[i]));
	}

	
	qsort(indices, n, sizeof(ItemIndex), compare_item_index_ab);

	
	double* new_a = (double*)malloc(n * sizeof(double));
	double* new_b = (double*)malloc(n * sizeof(double));
	double* new_p = (double*)malloc(n * sizeof(double));
	double* new_pw = (double*)malloc(n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		size_t idx = indices[i].index;
		new_a[i] = preInstances.a_ptr[idx];
		new_b[i] = preInstances.b_ptr[idx];
		new_p[i] = preInstances.p_ptr[idx];
		new_pw[i] = preInstances.p_weight[idx];
	}

	
	free(preInstances.a_ptr);
	free(preInstances.b_ptr);
	free(preInstances.p_ptr);
	free(preInstances.p_weight);
	preInstances.a_ptr = new_a;
	preInstances.b_ptr = new_b;
	preInstances.p_ptr = new_p;
	preInstances.p_weight = new_pw;
}
// instance data
void sort_instance_by_profit(Instance& instance) {
	size_t n = instance.n_items;
	ItemIndex* indices = (ItemIndex*)malloc(n * sizeof(ItemIndex));

	for (size_t i = 0; i < n; i++) {
		indices[i].index = i;
		indices[i].sum_p_ab = instance.p_ptr[i];
	}

	qsort(indices, n, sizeof(ItemIndex), compare_item_index_prt);

	double* new_a = (double*)malloc(n * sizeof(double));
	double* new_b = (double*)malloc(n * sizeof(double));
	double* new_p = (double*)malloc(n * sizeof(double));
	double* new_pw = (double*)malloc(n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		size_t idx = indices[i].index;
		new_a[i] = instance.a_ptr[idx];
		new_b[i] = instance.b_ptr[idx];
		new_p[i] = instance.p_ptr[idx];
		new_pw[i] = instance.p_weight[idx];
	}

	instance.a_ptr = new_a;
	instance.b_ptr = new_b;
	instance.p_ptr = new_p;
	instance.p_weight = new_pw;

	free(indices);
}

//convert vector to string
std::string JoinVector(const std::vector<int>& vec) {
	std::ostringstream oss;
	for (size_t i = 0; i < vec.size(); ++i) {
		if (i != 0) oss << ",";
		oss << vec[i];
	}
	return oss.str();
}

int FKP(const Instance* RESTRICT const inst,
	DblMatrix* RESTRICT const output,
	int binCapacity,
	const vector<double>& weightRec,
	bool dualFlag) {

	ItemList item_list = { 0 };
	run_or_fail(item_list_alloc(&item_list, inst, FullWeight));
	if (!dualFlag)
		invert_item_list(&item_list);

	size_t n_rows = output->cols;
	item_list.n_items = 0;
	size_t n_cols = weightRec.size() - 1;

	for (size_t i = 1; i < n_rows; ++i) {
		++item_list.n_items;
		item_list_sort(&item_list, AWeight);

		double last_w = 0;
		double last_p = 0;
		size_t start_item = 0;
		bool complete = false;
		for (size_t j = 0; j < n_cols; ++j) {
			double w = weightRec[j + 1];
			double curr_w = last_w;
			double curr_p = last_p;
			bool run = true;
			if (!complete) {
				for (size_t k = start_item; k < i && run; ++k) {
					Item item = item_list.item_ptr[k];
					double avail = w - curr_w;
					double inW = item.a_val;
					double p = item.p_val;
					if (avail < inW) {
						double frac = avail / inW;
						inW *= frac;
						p *= frac;
						run = false;
						start_item = k;
					}
					else {
						last_w += inW;
						last_p += p;
					}
					curr_w += inW;
					curr_p += p;
				}
				if (run)
					complete = true;
			}
			dm_set(output, i, j, curr_p);
		}
	}

	item_list_free(&item_list);
	//free(weight_to_secondIdx);
	return SUCCESS;
}

//label extention
void LabelExtention(
	Instance& preInstance,
	MyLabel* parentLab,
	MyLabel* preLab,
	int item
) {
	preLab->lastItem = item;
	preLab->sum_a += preInstance.a_ptr[item];
	preLab->sum_b += preInstance.b_ptr[item];
	preLab->tolWeight = preLab->sum_a + preInstance.rho * sqrt(preLab->sum_b);
	preLab->tolProfit += preInstance.p_ptr[item];
	preLab->parentLab = parentLab;
}

//Judge same label
bool JudgeEquality(MyLabel* lab1, MyLabel* lab2) {
	return abs(lab1->tolProfit - lab2->tolProfit) <= EX &&
		abs(lab1->tolWeight - lab2->tolWeight) <= EX &&
		abs(lab1->sum_a - lab2->sum_a) <= EX;
}

//dominance rules
bool JgeLabDominance(MyLabel* lab1, MyLabel* lab2, bool heuDom) {
	if (!heuDom) {//exact dominance rules
		if (!g_pInteger) {
			if (lab1->tolProfit < lab2->tolProfit - EX) return false;
		}
		else {
			if ((int)lab1->tolProfit < (int)lab2->tolProfit) return false;
		}
		//if (lab1->sum_b > lab2->sum_b + EX) return false; //old weaker version
		if (lab1->tolWeight > lab2->tolWeight + EX) return false;
		if (!g_aInteger) {
			if (lab1->sum_a > lab2->sum_a + EX) return false;
		}
		else {
			if ((int)lab1->sum_a > (int)lab2->sum_a) return false;
		}


		if (g_pInteger && !g_aInteger) {
			if (
				(int)lab1->tolProfit > (int)lab2->tolProfit ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				lab1->sum_a < lab2->sum_a - EX) return true;
		}
		else if (!g_pInteger && g_aInteger) {
			if (
				lab1->tolProfit > lab2->tolProfit + EX ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				(int)lab1->sum_a < (int)lab2->sum_a) return true;
		}
		else if (!g_pInteger && !g_aInteger) {
			if (
				lab1->tolProfit > lab2->tolProfit + EX ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				lab1->sum_a < lab2->sum_a - EX) return true;
		}
		else {//(g_pInteger && g_aInteger)
			if (
				(int)lab1->tolProfit > (int)lab2->tolProfit ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				(int)lab1->sum_a < (int)lab2->sum_a) return true;
		}

		return false;
	}
	else {//heurisitc dominance rules
		if (lab1->tolProfit < lab2->tolProfit - EX) return false;
		if (lab1->tolWeight > lab2->tolWeight + EX) return false;

		if (
			lab1->tolProfit > lab2->tolProfit + EX ||
			lab1->tolWeight < lab2->tolWeight - EX
			) return true;
		return false;
	}
}

//completion bound
int CompletionBound(
	MyLabel* lab,
	int prefirstIdx,
	DblMatrix* ub_matr
) {
	//strengthened completion bound
	size_t i = ub_matr->cols - prefirstIdx - 1;
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];
	size_t j = ub_matr->rows - secondIdx - 1;

	if (i < g_firstDim) {
		double leftVal = dm_get(ub_matr, i, j) + lab->tolProfit;
		if (leftVal <= g_bestLB + EX) {//
			return true;
		}
	}
	return false;
}

//label dominance
void JgeDominance(
	vector<vector<multimap<double, MyLabel*, greater<double>>>>& nonExtended,
	//vector<vector<multimap<double, MyLabel*, greater<double>>>>& extended,
	set<int>& firstNonBucketIdx,
	vector<set<int>>& secondNonBucketIdx,
	MyLabel* lab,
	int& labelCnt) {
	//find the dimension of the bucket
	int firstIdx = g_item_to_firstIdx[lab->lastItem];
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];

	////judge dominance
	bool dominanceFlag = false;
	auto idxIte1 = firstNonBucketIdx.begin();

	if (!dominanceFlag) {
		//use the present label to dominate other labels in the bucket
		idxIte1 = firstNonBucketIdx.find(firstIdx);
		while (idxIte1 != firstNonBucketIdx.end()) {
			int i = *idxIte1;
			auto idxIte2 = secondNonBucketIdx[i].find(secondIdx);
			auto tIte = idxIte2;
			while (idxIte2 != secondNonBucketIdx[i].end()) {
				int j = *idxIte2;
				if (!nonExtended[i][j].empty() &&
					lab->tolProfit < (--nonExtended[i][j].end())->first) {
					++idxIte2;
					continue;
				}

				auto domIite = nonExtended[i][j].begin();
				while (domIite != nonExtended[i][j].end()) {
					if (JgeLabDominance(lab, domIite->second, true)) {
						delete domIite->second;
						domIite = nonExtended[i][j].erase(domIite);
						--labelCnt;
						//cout << "nonExtended dominanted" << endl;
						continue;
					}
					++domIite;
				}
				
				++idxIte2;
			}
			++idxIte1;
		}
	}
}

//label dominance
void JgeDominance(
	vector<multimap<double, MyLabel*, greater<double>>>& newExtended,
	vector<multimap<double, MyLabel*, greater<double>>>& oldExtended,
	MyLabel* lab,
	vector<MyLabel*>& dominatedOldLabs,
	bool newLabelFlag = true,
	bool heuDom = false,
	bool dualCalFlag = false) {
	//find the dimension of the bucket
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];

	//use the present label to dominate other labels in the bucket
	int j = secondIdx;
	while (j < secondDim) {
		//dominate the labels in the new bucket
		if (!newExtended[j].empty() &&
			lab->tolProfit >= (--newExtended[j].end())->first) {
			auto ite = newExtended[j].begin();
			while (ite != newExtended[j].end()) {
				if (JgeLabDominance(lab, ite->second, heuDom)) {
					if (dualCalFlag) {
						delete ite->second;
						ite->second = nullptr;
					}
					else {
						if (j == secondIdx)
							dominatedOldLabs.push_back(ite->second);
						else {
							delete ite->second;
							ite->second = nullptr;
							++g_dominatedLabel;
						}
					}
					ite = newExtended[j].erase(ite);
					continue;
				}
				else {//Judge label equality
					if (JudgeEquality(lab, ite->second)) {
						if (dualCalFlag) {
							delete ite->second;
							ite->second = nullptr;
						}
						else {
							if (j == secondIdx)
								dominatedOldLabs.push_back(ite->second);
							else {
								delete ite->second;
								ite->second = nullptr;
								++g_dominatedLabel;
							}
						}
						ite = newExtended[j].erase(ite);
						continue;
					}
				}
				++ite;
			}
		}

		++j;
		if (j - secondIdx > g_domWidth)
			break;
	}
}


//get the initial lower bound
double GetInitialLowerBound(double& bestLB) {
	Instance instance = g_instance;
	//instance profit
	sort_instance_by_profit(instance);

	double sum_a = 0;
	double sum_b = 0;
	for (int i = 0; i < instance.n_items; ++i) {
		sum_a += instance.a_ptr[i];
		sum_b += instance.b_ptr[i];
		if (sum_a + instance.rho * sqrt(sum_b) > instance.capacity)
			break;
		else
			bestLB += instance.p_ptr[i];
	}

	MyLabel* preLab = new MyLabel();
	preLab->tolProfit = bestLB;
	return bestLB;
}

//initialize g_profitRec
void InitProfitRec() {
	vector<int> tmp(g_instance.n_items, 0);
	g_profitRec.resize(g_instance.n_items, 0);
	double tolProfit = 0;
	tmp[0] = g_instance.p_ptr[0];
	for (int i = 0; i < g_instance.n_items; ++i) {
		tolProfit += g_instance.p_ptr[i];
		if (i >= 1)
			tmp[i] = tmp[i - 1] + g_instance.p_ptr[i];
	}
	for (int i = 0; i < g_instance.n_items; ++i)
		g_profitRec[i] = tolProfit - tmp[i];
}

//The whole dominance logic
bool DominanceLogic(
	int preSecondIdx,
	vector<multimap<double, MyLabel*, greater<double>>>*& newExtended,
	vector<multimap<double, MyLabel*, greater<double>>>*& oldExtended,
	MyLabel* preLab,
	bool newLabelFalg = true,
	bool heuDom = false) {
	bool dominanceFlag = false;
	int j = preSecondIdx;
	while (j >= 0) {
		auto& preBkt = (*newExtended)[j];
		//use newExtended labels to dominante
		auto domIite = preBkt.begin();
		int cnt = 0;
		int preSize = preBkt.size();
		while (cnt < preSize) {
			if (domIite->first < preLab->tolProfit) break;
			//judge dominance
			if (JgeLabDominance(domIite->second, preLab, heuDom)) {
				dominanceFlag = true;
				break;
			}
			else {
				//judge equality
				if (JudgeEquality(domIite->second, preLab)) {
					dominanceFlag = true;
					//cout << "Equal label 1" << endl;
					break;
				}
			}
			++domIite;
			++cnt;
		}
		if (dominanceFlag) break;

		//old labels don't need to be checked from old labels, since it's checked during the geration of the new label
		if (newLabelFalg) {
			auto& preBkt = (*oldExtended)[j];
			auto domIite = preBkt.begin();
			int cnt = 0;
			int preSize = preBkt.size();
			//while (domIite != preBkt.end()) {
			while (cnt < preSize) {
				if (domIite->first < preLab->tolProfit) break;
				//judge dominance
				if (JgeLabDominance(domIite->second, preLab, heuDom)) {
					dominanceFlag = true;
					break;
				}
				else {
					//judge equality
					if (JudgeEquality(domIite->second, preLab)) {
						dominanceFlag = true;
						//cout << "Equal label 2" << endl;
						break;
					}
				}
				++domIite;
				++cnt;
			}
			if (dominanceFlag)break;
		}

		--j;
		if (preSecondIdx - j > g_domWidth) break;
	}

	return dominanceFlag;
}


//use the heuristic labeling algorithm with heuristic dominance rules to solve the submodular knapsack problem
void LabelSettingHeuristic(
	DblMatrix& ub_matr,
	double& heuLableTime,
	const ItemIndex* indicesRec,
	vector<double>& SR3Duals,
	vector<unordered_set<int>>& SR3s,
	unordered_map<int, ItemIndex>& removedIndicesRec,
	multimap<double, KnapsackSol, greater<double>>& finalSols
	//ofstream& outPut
) {
	auto startTime = chrono::high_resolution_clock::now();

	//initialize the bucket
	vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
	dominatedOldLabs.reserve(1e+6);
	vector<int> nonDominatedLabsRec(g_firstDim);

	MyLabel* initLab = new MyLabel();
	oldExtended->begin()->insert({ 0, initLab });
	++g_tolGeneratedLabel;

	//label extention and dominance
	int currItem = 0;				//record the current item to extend
	int thisNonDominatedLabel = 0;
	for (int stage = 0; stage < g_firstDim; ++stage) {
		thisNonDominatedLabel = 0;
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*oldExtended)[sndBktIdx].empty()) continue;
			auto ite = (*oldExtended)[sndBktIdx].begin();
			while (ite != (*oldExtended)[sndBktIdx].end()) {
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
					dominatedOldLabs.push_back(ite->second);
					keepOld = false;
				}
				else {
					++g_tolGeneratedLabel;
				}

				double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
					g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
				if (preWeight <= g_instance.capacity + EX) {//judge capacity
					//label extension
					MyLabel* tmpLab = new MyLabel(ite->second);
					LabelExtention(g_instance, ite->second, tmpLab, currItem);
					++g_tolGeneratedLabel;

					//completion bound to fathom label
					if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
						++g_CBFathomLabel;
						delete tmpLab; tmpLab = nullptr;
						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, true);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
								continue;
							}
							(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
						}
						continue;
					}

					//dominance check
					int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
					bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab, true, true);

					if (dominanceFlag) {
						delete tmpLab;
						tmpLab = nullptr;
						++g_dominatedLabel;
					}
					else {
						//dominate other labels
						JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs, true, true);
						//add the non dominated label into the new extended bucket
						(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*oldExtended)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, true);
					if (dominanceFlag) {
						/*delete oldLab;
						oldLab = nullptr;*/
						dominatedOldLabs.push_back(oldLab);
						continue;
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
			}
			thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
			//use the present best dual variable to update
			if (!g_solveKnapsack && !(*newExtended)[sndBktIdx].empty()) {
				if ((*newExtended)[sndBktIdx].begin()->first > g_bestLB) {
					KnapsackSol preSol;
					preSol.bestLab = new MyLabel((*newExtended)[sndBktIdx].begin()->second);
					//get the item set
					if (preSol.bestLab->itemSet.empty())
						preSol.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
					auto tmpLab = preSol.bestLab;
					while (tmpLab->parentLab != nullptr) {
						preSol.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
						tmpLab = tmpLab->parentLab;
					}
					finalSols.insert({ (*newExtended)[sndBktIdx].begin()->first, preSol });
					g_bestLB = max(g_bestLB, (*newExtended)[sndBktIdx].begin()->first);
				}
			}
		}

		++currItem;
		delete oldExtended;
		oldExtended = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		nonDominatedLabsRec[stage] = thisNonDominatedLabel;
	}

	//record the best solution
	auto ite = --(*oldExtended).end();
	MyLabel* bestLab = nullptr;
	while (true) {
		if (!ite->empty()) {
			if (bestLab == nullptr)
				bestLab = ite->begin()->second;

			auto subIte = ite->begin();
			while (subIte != ite->end()) {
				/*if (1 - subIte->first >= 0)
					break;*/
				KnapsackSol preSol;
				preSol.bestLab = subIte->second;
				finalSols.insert({ subIte->second->tolProfit, preSol });
				subIte = ite->erase(subIte);
			}
		}
		if (ite == (*oldExtended).begin())
			break;
		--ite;
	}
	if (bestLab != nullptr && bestLab->lastItem > 0) {
		//get the item set
		vector<int> bestIS;
		MyLabel* tmpLab = bestLab;
		while (tmpLab->parentLab != nullptr) {
			bestIS.push_back(indicesRec[tmpLab->lastItem].index);
			tmpLab = tmpLab->parentLab;
		}
		sort(bestIS.begin(), bestIS.end());
		string bestItemSet = JoinVector(bestIS);
		finalSols.begin()->second.bestItemSet = bestItemSet;
		//get the item set for all the solutions
		for (auto& e : finalSols) {
			if (e.second.bestLab->itemSet.empty())
				e.second.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
			tmpLab = e.second.bestLab;
			while (tmpLab->parentLab != nullptr) {
				e.second.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
				tmpLab = tmpLab->parentLab;
			}
		}
	}
	//free space
	for (auto& t : (*oldExtended)) {
		for (auto& e : t)
			delete e.second;
	}
	delete oldExtended;
	delete newExtended;
	g_dominatedLabel += dominatedOldLabs.size();
	for (auto& t : dominatedOldLabs)
		delete t;
	g_nonDominatedLabel = thisNonDominatedLabel;
	//cout << "The trend of non dominated labels in heuristic label algorithm: " << endl;
	//outPut << "The trend of non dominated labels in heuristic label algorithm: " << endl;
	//for (auto& e : nonDominatedLabsRec) {
	//	cout << e << "\t";
	//	outPut << e << "\t";
	//}
	cout << endl;
	//outPut << endl;
	auto endTime = chrono::high_resolution_clock::now();
	heuLableTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
}

//test solution
bool TestSolution(Instance& preInstance, vector<int>& sol_items) {
	double sum_a = 0;
	double sum_b = 0;
	double profit = 0;
	for (auto& e : sol_items) {
		if (e > 0) {
			sum_a += preInstance.a_ptr[e];
			sum_b += preInstance.b_ptr[e];
			profit += preInstance.p_ptr[e];
		}
	}
	if (sum_a + preInstance.rho * sqrt(sum_b) <= preInstance.capacity) {
		cout << "Solution feasible, with profit: " << profit << endl;
		return true;
	}
	else {
		cout << "Solution infeasible! " << endl;
		return false;
	}
}

//design the labelsetting algorithm to solve the submodular knapsack problem
int LabelSettingSolveKnapsack(Args& args) {
	auto startTime = chrono::high_resolution_clock::now();
	/*read g_instance*/
	g_instance = { 0 };
	if (instance_parse(&g_instance, args.input_file)) {
		cerr << "Error reading g_instance file" << endl;
		return 0;
	}
	ofstream outPut(args.output_file, ios::app);
	//judge if p and a is integer or not
	for (int i = 0; i < g_instance.n_items; ++i) {
		if (abs(floor(g_instance.a_ptr[i]) - g_instance.a_ptr[i]) > EX) {
			g_aInteger = false;
		}
		if (abs(floor(g_instance.p_ptr[i]) - g_instance.p_ptr[i]) > EX) {
			g_pInteger = false;
		}
		double preRatio = g_instance.b_ptr[i] * 1.0 / (g_instance.a_ptr[i] * 1.0);
		if (preRatio < g_smallest_bDiva)
			g_smallest_bDiva = preRatio;
	}

	//specify the dimension of the bucket
	g_tolItemNum = g_instance.n_items;			//the first dimension of the bucket
	g_firstDim = g_instance.n_items / 1;		//the first dimension of the bucket
	secondDim = g_instance.capacity / g_secondDimDiv;       //the second dimension of the bucket

	/*use heuristic to get primal bound*/
	auto heuPrimal_startTime = chrono::high_resolution_clock::now();
	//get the initial lower bound
	multimap<double, KnapsackSol, greater<double>> finalSols;
	g_bestLB = 0;
	GetInitialLowerBound(g_bestLB);
	//sort the instance according to p/(a+sqrt(b))
	ItemIndex* indicesRec = nullptr;
	sort_instance_by_p_ab(g_instance, indicesRec);
	//TS heuristic
	vector<int> ts_sol_items(g_tolItemNum);
	double ts_sol = heuristic_solution(&g_instance, &ts_sol_items[0]);
	g_bestLB = max(g_bestLB, ts_sol);
	double tsSol_sum_a = 0;
	double tsSol_sum_b = 0;
	vector<int> map_ts_sol_items(g_tolItemNum);
	for (int i = 0; i < g_tolItemNum; ++i) {
		if (ts_sol_items[i]) {
			tsSol_sum_a += g_instance.a_ptr[i];
			tsSol_sum_b += g_instance.b_ptr[i];
			map_ts_sol_items[indicesRec[i].index] = 1;
		}
	}
	MyLabel* tsSolLab = new MyLabel();
	tsSolLab->tolProfit = ts_sol;
	tsSolLab->itemSet = map_ts_sol_items;
	tsSolLab->tolWeight = tsSol_sum_a + g_instance.rho * sqrt(tsSol_sum_b);
	KnapsackSol preSol;
	preSol.bestLab = tsSolLab;
	finalSols.insert({ ts_sol, preSol });
	auto heuPrimal_endTime = chrono::high_resolution_clock::now();

	/*construct bucket index*/
	int itemGap = g_instance.n_items / g_firstDim;	//the item gap for each item bucket
	g_item_to_firstIdx.resize(g_instance.n_items + 1, 0);
	for (int j = 0; j <= g_firstDim; ++j) {
		int startItem = j * itemGap;
		int endItem = (j == g_firstDim - 1) ? g_instance.n_items : (j + 1) * itemGap - 1;
		endItem = min((size_t)endItem, g_instance.n_items);
		for (int w = startItem; w <= endItem; ++w)
			g_item_to_firstIdx[w] = j;
	}
	vector<double> weightRec(secondDim + 1, 0);
	int weightGap = g_instance.capacity / secondDim;	//the capacity gap for each capacity bucket
	g_weight_to_secondIdx.resize(g_instance.capacity + 1, 0);
	for (int j = 0; j <= secondDim; ++j) {
		int startWeight = j * weightGap;
		int endWeight = (j == secondDim - 1) ? g_instance.capacity : (j + 1) * weightGap - 1;
		endWeight = min(endWeight * 1.0, g_instance.capacity);
		weightRec[j] = endWeight;
		for (int w = startWeight; w <= endWeight; ++w)
			g_weight_to_secondIdx[w] = j;
	}

	/*compute linear upper bound (dual bound)*/
	auto LR_startTime = chrono::high_resolution_clock::now();
	DblMatrix ub_matr = { 0 };
	dm_new(&ub_matr, g_instance.n_items + 1, secondDim + 1);
	for (size_t i = 0; i <= g_instance.n_items; ++i) {
		for (size_t j = 0; j <= secondDim; ++j) {
			dm_set(&ub_matr, i, j, INFINITY);
		}
	}
	FKP(&g_instance, &ub_matr, g_instance.capacity + 1, weightRec, false);//mask completion bound
	auto LR_endTime = chrono::high_resolution_clock::now();

	/*use the heuristic labeling algorithm to get better lower bound*/
	double heuLableTime = 0;
	unordered_map<int, ItemIndex> removedIndicesRec;
	vector<double> sr3Duals;
	vector<unordered_set<int>> sr3s;
	vector<double> preciseDuals(g_instance.n_items, 0);
	for (int i = 0; i < g_instance.n_items; ++i)
		preciseDuals[i] = g_instance.p_ptr[i];

	double dualBoundTime = 0; bool solvedFlag = true;

	LabelSettingHeuristic(ub_matr, heuLableTime, indicesRec, sr3Duals, sr3s, removedIndicesRec, finalSols);//mask HLA
	cout << "heuristic labeling get lower bound: " << g_bestLB << endl;
	cout << "heuristic labeling use time: " << heuLableTime << "s" << endl;
	int HLA_sol = g_bestLB;
	
	auto exactDP_startTime = chrono::high_resolution_clock::now();
	//initialize the bucket
	vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
	dominatedOldLabs.reserve(1e+6);
	vector<int> nonDominatedLabsRec(g_firstDim);

	MyLabel* initLab = new MyLabel();
	oldExtended->begin()->insert({ 0, initLab });
	++g_tolGeneratedLabel;

	//label extention and dominance
	int currItem = 0;				//record the current item to extend
	bool timeLimitFlag = false;
	int thisNonDominatedLabel = 0;
	for (int stage = 0; stage < g_firstDim; ++stage) {
		thisNonDominatedLabel = 0;
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*oldExtended)[sndBktIdx].empty()) continue;
			auto ite = (*oldExtended)[sndBktIdx].begin();
			while (ite != (*oldExtended)[sndBktIdx].end()) {
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
					dominatedOldLabs.push_back(ite->second);
					keepOld = false;
				}
				else {
					++g_tolGeneratedLabel;
				}

				double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
					g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
				if (preWeight <= g_instance.capacity + EX) {//judge capacity
					//label extension
					MyLabel* tmpLab = new MyLabel(ite->second);
					LabelExtention(g_instance, ite->second, tmpLab, currItem);
					++g_tolGeneratedLabel;

					//completion bound to fathom label
					if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
						++g_CBFathomLabel;
						delete tmpLab; tmpLab = nullptr;
						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
								continue;
							}

							(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
						}
						continue;
					}

					//dominance check
					int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
					bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab);

					if (dominanceFlag) {
						delete tmpLab;
						++g_dominatedLabel;
					}
					else {
						//dominate other labels
						JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs);
						//add the non dominated label into the new extended bucket
						(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*oldExtended)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
					if (dominanceFlag) {
						/*delete oldLab;
						oldLab = nullptr;*/
						dominatedOldLabs.push_back(oldLab);
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
			}
			thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
		}
		auto endTime = chrono::high_resolution_clock::now();
		if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
			cout << "The maximum solution time of my DP is reached!" << endl;
			outPut << "The maximum solution time of my DP is reached!" << endl;
			timeLimitFlag = true;
		}
		if (timeLimitFlag) {
			auto ite = (*oldExtended).begin();
			for (auto& e : (*newExtended)) {
				for (auto& t : e)
					ite->insert({ t.first, t.second });
				++ite;
			}
			break;
		}
		++currItem;
		delete oldExtended;
		oldExtended = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		nonDominatedLabsRec[stage] = thisNonDominatedLabel;
	}
	//record the best solution org record for ELA
	auto ite = --(*oldExtended).end();
	MyLabel* bestLab = nullptr;
	while (true) {
		if (!ite->empty()) {
			bestLab = ite->begin()->second;
			ite->erase(ite->begin());
			break;
		}
		if (ite == (*oldExtended).begin())
			break;
		--ite;
	}
	string bestItemSet = finalSols.begin()->second.bestItemSet;
	if (bestLab != nullptr) {
		//get the item set
		vector<int> bestIS;
		MyLabel* tmpLab = bestLab;
		while (tmpLab->parentLab != nullptr) {
			bestIS.push_back(indicesRec[tmpLab->lastItem].index);
			tmpLab = tmpLab->parentLab;
		}
		sort(bestIS.begin(), bestIS.end());
		bestItemSet = JoinVector(bestIS);
		//store the best solution
		KnapsackSol bestSol;
		bestSol.bestItemSet = bestItemSet;
		bestSol.bestLab = bestLab;
		finalSols.insert({ bestLab->tolProfit, bestSol });
		g_bestLB = max(g_bestLB, finalSols.begin()->first);
	}
	if (bestItemSet.empty()) {
		for (int i = 0; i < finalSols.begin()->second.bestLab->itemSet.size(); ++i)
			if (finalSols.begin()->second.bestLab->itemSet[i])
				bestItemSet += to_string(i) + ",";
	}

	//free space
	for (auto& e : finalSols)
		delete e.second.bestLab;
	for (auto& t : (*oldExtended)) {
		for (auto& e : t)
			delete e.second;
	}
	delete oldExtended;
	delete newExtended;
	g_dominatedLabel += dominatedOldLabs.size();
	for (auto& lab : dominatedOldLabs)
		delete lab;

	free(indicesRec);
	free(g_instance.a_ptr);
	free(g_instance.b_ptr);
	free(g_instance.p_ptr);
	free(g_instance.p_weight);
	g_instance.a_ptr = g_instance.b_ptr = g_instance.p_ptr = g_instance.p_weight = nullptr;
	dm_free(&ub_matr);

	auto exactDP_endTime = chrono::high_resolution_clock::now();
	auto endTime = chrono::high_resolution_clock::now();
	double myLableTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
	cout << "The best solution obtained by my labelsetting algorithm is: " << g_bestLB << endl;
	cout << "The solution obtained by TS heuristic is: " << ts_sol << endl;
	cout << "The solution obtained by heuristic labeling algorithm is: " << HLA_sol << endl;
	cout << "The time for my labelsetting algorithm is: " << myLableTime << "s" << endl;
	cout << "The time for heuristic primal solution is: " << chrono::duration_cast<chrono::milliseconds>(heuPrimal_endTime - heuPrimal_startTime).count() / 1000.0 << "s" << endl;
	cout << "The time for linear relaxation (dual bound) is: " << chrono::duration_cast<chrono::milliseconds>(LR_endTime - LR_startTime).count() / 1000.0 << "s" << endl;
	cout << "The time for heuristic DP is: " << heuLableTime << "s" << endl;
	cout << "The time for exact dual bound is: " << dualBoundTime << "s" << endl;
	cout << "The time for final exact DP is: " << chrono::duration_cast<chrono::milliseconds>(exactDP_endTime - exactDP_startTime).count() / 1000.0 << "s" << endl;
	cout << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
	cout << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
	cout << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
	cout << "The number of dominated labels is: " << g_dominatedLabel << endl;
	cout << "The best item set is: " << bestItemSet << endl;
	//cout << "The trend of non dominated labels: " << endl;

	outPut << "\n\n****The resulte of my labelsetting algorithm****" << g_bestLB << endl;
	outPut << "The best solution obtained by my labelsetting algorithm is: " << g_bestLB << endl;
	outPut << "The solution obtained by TS heuristic is: " << ts_sol << endl;
	outPut << "The solution obtained by heuristic labeling algorithm is: " << HLA_sol << endl;
	outPut << "The time for my labelsetting algorithm is: " << myLableTime << "s" << endl;
	outPut << "The time for heuristic primal solution is: " << chrono::duration_cast<chrono::milliseconds>(heuPrimal_endTime - heuPrimal_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The time for linear relaxation (dual bound) is: " << chrono::duration_cast<chrono::milliseconds>(LR_endTime - LR_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The time for heuristic DP is: " << heuLableTime << "s" << endl;
	outPut << "The time for exact dual bound is: " << dualBoundTime << "s" << endl;
	outPut << "The time for final exact DP is: " << chrono::duration_cast<chrono::milliseconds>(exactDP_endTime - exactDP_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
	outPut << "The best item set is: " << bestItemSet << endl;
	outPut << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
	outPut << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
	outPut << "The number of dominated labels is: " << g_dominatedLabel << endl;
	//outPut << "The trend of non dominated labels: " << endl;
	//for (auto& e : nonDominatedLabsRec) {
	//	cout << e << "\t";
	//	outPut << e << "\t";
	//}
	cout << endl;
	outPut << endl;
	outPut.close();

	return g_bestLB;
}