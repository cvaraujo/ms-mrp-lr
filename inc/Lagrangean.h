//
// Created by carlos on 31/05/19.
//

#ifndef LAGRANGEANMST_LAGRANGEAN_H
#define LAGRANGEANMST_LAGRANGEAN_H


#include "Graph.h"
#include "Model.h"
#include <chrono>

class Lagrangean {
    Graph *graph;
    Model *model;
    double thetaC, thetaP, thetaD, objectiveFunction, originalObjective, lambda, LB;
    int progress, iter, maxIter, B, UB, time, iterBlb, iterBub, endTime;
    bool feasible = true;
    vector<double> multipliersDelay, multipliersJitter;
    vector<vector<double >> multipliersVar;
    vector<vector<vector<double >>> multipliersRel;

    bool solveModel();

public:
    Lagrangean(Graph *graph, double lambda, int maxIter, int B, int time);

    void getGradient(vector<double> &gradientDelay, vector<double> &gradientJitter);

    void getGradientTerminals(vector<vector<double>> &gradientVar);

    void getGradientRelation(vector<vector<vector<double>>> &gradientRel);

    double getNormTerminals(vector<vector<double>> &gradient);

    double getNormTerminals(vector<double> &gradient);

    double getNormRelation(vector<vector<vector<double>>> &gradient);

    void updatePPL();

    void updatePathCostsNT();

    double originalObjectiveValue();

    bool isFeasible();

    void updateTreeCosts();

    void updatePathCosts(int k);

    int heuristic();

    double solve();

    void showSolution(string outputName);
};


#endif //LAGRANGEANMST_LAGRANGEAN_H
