//
// Created by carlos on 31/05/19.
//

#ifndef LAGRANGEANMST_LAGRANGEAN_H
#define LAGRANGEANMST_LAGRANGEAN_H


#include "Graph.h"
#include "Model.h"
#include "BarrierMethod.h"
#include <chrono>

class Lagrangean {
  Graph *graph;
  Model *model;
  double thetaC, thetaP, thetaD, objectiveFunction, originalObjective, lambda, LB, firstLB;
  int progress, iter, maxIter, B, UB, time, iterBlb, iterBub, bmTime = 0, endTime, relaxNum, firstUB;
  
  bool feasible = true, heuristics = false;
  vector<double> multipliersDelay, multipliersJitter;
  vector<vector<double >> multipliersVar, multipliersLeaf;
  vector<vector<vector<double >>> multipliersRel;
  vector<vector<int>> freqLB, freqUB;

  bool solveModel();

  void getGradientDelay(vector<double> &gradientDelay);

  void getGradientJitter(vector<double> &gradientJitter);

  void getGradientVariation(vector<vector<double>> &gradientVar);

  void getGradientLeaf(vector<vector<double>> &gradientLeaf);

  void getGradientRelation(vector<vector<vector<double>>> &gradientRel);

  double getNormDelay(vector<double> &gradientDelay);

  double getNormJitter(vector<double> &gradientJitter);

  double getNormRelation(vector<vector<vector<double>>> &gradient);

  double getNormLeaf(vector<vector<double>> &gradientLeaf);

  double getNormVariation(vector<vector<double>> &gradient);

public:
  Lagrangean(Graph *graph, int relaxNum, bool heuristics, bool barrierMethod, double lambda, int maxIter, int B, int time);



  double getNormTerminals(vector<double> &gradient);



  void updatePPL();

  void updatePathCostsNT(int q);

  double originalObjectiveValue();

  bool isFeasible();

  void updateTreeCosts();

  void updatePathCosts(int k);

  int heuristic();

  double solve();

  void showSolution(string outputName);
};


#endif //LAGRANGEANMST_LAGRANGEAN_H
