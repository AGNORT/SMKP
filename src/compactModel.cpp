#include "compactModel.h"
#include "labelSetting.h"

using namespace std;

//call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
void SolveCompactKnapsackModel(Args& args) {
    //read instance
    Instance instance = { 0 };
    if (instance_parse(&instance, args.input_file)) {
        cerr << "Error reading instance file" << endl;
        return;
    }
    ofstream outPut(args.output_file, ios::app);

    try {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "./gurobi.log");
        env.start();
        GRBModel model = GRBModel(env);

        //variable specify if the item i is selected
        GRBVar* x = model.addVars(instance.n_items, GRB_BINARY);
        GRBVar t = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        //set objective function
        GRBLinExpr obj = 0;
        for (int i = 0; i < instance.n_items; ++i)
            obj += instance.p_ptr[i] * x[i];
        model.setObjective(obj, GRB_MAXIMIZE);  

        //set SOCP constraints
        GRBLinExpr lhs = 0;
        for (int j = 0; j < instance.n_items; ++j)
            lhs += instance.a_ptr[j] * x[j];
        model.addConstr(lhs + instance.rho * t <= instance.capacity, "SOCP_constraint");
        GRBLinExpr rhs = 0;
        for (int j = 0; j < instance.n_items; ++j)
            rhs += instance.b_ptr[j] * x[j];
        model.addQConstr(t * t == rhs, "t_squared_constraint");

        // set precision
        model.set(GRB_DoubleParam_MIPGap, 1e-6);           
        model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);   
        model.set(GRB_DoubleParam_OptimalityTol, 1e-9);    
        model.set(GRB_DoubleParam_IntFeasTol, 1e-9);       
        model.set(GRB_IntParam_Presolve, 2);               
        model.set(GRB_IntParam_Method, 2);                 

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);        //set maximum solution time

        // solve the model
        model.optimize();

        cout << "***************The results obtained by Gurobi**************" << endl;
        outPut << "\n\n***************The results obtained by Gurobi**************" << endl;
        // output the results
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout << "The optimal objective value is :" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            cout << "The solution time is :" << model.get(GRB_DoubleAttr_Runtime) << endl;
            cout << "The chosen items are as follows: " << endl;
            outPut << "The optimal objective value is :" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            outPut << "The solution time is :" << model.get(GRB_DoubleAttr_Runtime) << endl;
            outPut << "The chosen items are as follows: " << endl;

            int itemNum = 0;
            for (int i = 0; i < instance.n_items; i++) {
                if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
                    ++itemNum;
                    cout << i << ",";
                    outPut << i << ",";
                }
            }
            cout << endl;
            cout << "The total number of chosen items is: " << itemNum << endl;
            outPut << endl;
            outPut << "The total number of chosen items is: " << itemNum << endl;
            std::cout << "t = " << t.get(GRB_DoubleAttr_X) << std::endl;
        }
        else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            cout << "Time limit reached. Returning the best solution found so far." << endl;
            outPut << "Time limit reached. Returning the best solution found so far." << endl;

            // Retrieve the best solution found up to the timeout
            double bestObjVal = model.get(GRB_DoubleAttr_ObjVal);
            cout << "The best objective value found is: " << bestObjVal << endl;
            outPut << "The best objective value found is: " << bestObjVal << endl;

            int itemNum = 0;
            for (int i = 0; i < instance.n_items; i++) {
                if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
                    ++itemNum;
                    cout << i << ",";
                    outPut << i << ",";
                }
            }
            cout << endl;
            cout << "The total number of chosen items is: " << itemNum << endl;
            outPut << endl;
            outPut << "The total number of chosen items is: " << itemNum << endl;
            std::cout << "t = " << t.get(GRB_DoubleAttr_X) << std::endl;
        }
        else {
            std::cout << "No optimal solution found" << std::endl;
            outPut << "No optimal solution found" << std::endl;
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi error code: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        outPut << "Gurobi error code: " << e.getErrorCode() << std::endl;
        outPut << e.getMessage() << std::endl;
        return;
    }
    catch (exception& e) {
        std::cout << "Other errors:" << e.what() << std::endl;
        outPut << "Other errors:" << e.what() << std::endl;
        return;
    }

    outPut.close();
    //free space
    instance_free(&instance);
}
