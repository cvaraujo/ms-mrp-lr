//
// Created by carlos on 06/03/19.
//

#ifndef MRP_BARRIER_H
#define MRP_BARRIER_H

#include "Graph.h"
#include "boost/algorithm/string.hpp"
#include <gurobi_c++.h>

class BarrierMethod {
    Graph *graph;
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    vector<vector<vector<GRBVar>>> f;
    vector<vector<GRBVar>> y;
    vector<GRBVar> z;
    int beginConstrRel = 0, endConstrRel = 0,
        beginConstrVar = 0, endConstrVar = 0,
        beginConstrLeaf = 0, endConstrLeaf = 0;


    void objectiveFunction();

    void rootFlow();

    void flowConservation();

    void terminalsFlow();

    void relFandY();

    void maxArcs();

    void limDelayAndJitter();

    void limVariation();

    void primeToTerminals();

    void nonTerminalsLeafs();

public:
	vector<double> multipliersDelay, multipliersJitter;
    vector<vector<double>> multipliersVar, multipliersLeaf;
    vector<vector<vector<double>>> multipliersRel;

    BarrierMethod(Graph *graph);

    void initialize();

    void initializeLinear();

    void initModel();
    
    void initModelLinearRelaxation();

    void initModelCshp();

    void solve();

    void solveLinear();

    void barrierMethod();

    void showSolution(string instance);

    int lagrangean();

};


#endif //MRP_MODEL_H
