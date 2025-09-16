#include<iostream>
#include "compactModel.h"
#include "labelSetting.h"
#include <string.h>

using namespace std;


int main(int argc, char** argv) {
	/**********Read arguments*******/
	Args args = { 0 };
	parse_args(&args, argc, argv);

	/**********Solve SMKP*******/
	cout << "Solve SMKP use ";
	if (strcmp(args.method, "DP") == 0) {
		cout << "DP!" << endl;
		//call labelsetting algorithm to solve the submodular knapsack problem
		LabelSettingSolveKnapsack(args);
	}
	else{
		cout << "Gurobi!" << endl;
		//call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
		SolveCompactKnapsackModel(args);
	}

	return 0;
}
