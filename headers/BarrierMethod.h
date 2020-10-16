//
// Created by carlos on 06/03/19.
//

#ifndef MRP_BARRIER_H
#define MRP_BARRIER_H

#include "Graph.h"
#include <gurobi_c++.h>

class BarrierMethod {
  Graph *graph;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  vector<vector<vector<GRBVar>>> f;
  vector<vector<GRBVar>> y;
  vector<GRBVar> z;

  void objectiveFunction();

  void rootFlow();

  void flowConservation();

  void terminalsFlow();

  void relXandY();

  void maxArcs();

  void limDelayAndJitter();

  void limVariation();

  void primeToTerminals();

  void nonTerminalsLeafs();

public:
  
  BarrierMethod(Graph *graph);

  void initialize();

  void initModel();
    
  void solve();

  void getMultipliersDelay(vector<double> &multipliersDelay);

  void getMultipliersJitter(vector<double> &multipliersJitter);

  void getMultipliersVariation(vector<vector<double>> &multipliersVar);

  void getMultipliersRelation(vector<vector<vector<double>>> &multipliersRel);

  void getMultipliersLeaf(vector<vector<double>> &multipliersLeaf);

  void showSolution(string instance);


};


#endif //MRP_MODEL_H
