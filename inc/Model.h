//
// Created by carlos on 28/05/19.
//

#ifndef LAGRANGEANMST_MODEL_H
#define LAGRANGEANMST_MODEL_H

#include "Graph.h"
#include "boost/graph/strong_components.hpp"

using namespace std;
using namespace boost;

struct SPPRC_Graph_Vert {
    SPPRC_Graph_Vert(int n = 0, int c_1 = 0, int c_2 = 0) : num(n), con_1(c_1), con_2(c_2) {}

    int num;
    // Resource consumed
    int con_1;
    int con_2;
};

struct SPPRC_Graph_Arc {
    SPPRC_Graph_Arc(int n = 0, double c = 0, int r_1 = 0, int r_2 = 0) : num(n), cost(c), res_1(r_1), res_2(r_2) {}

    int num;
    // traversal cost
    double cost;
    // traversal resource
    int res_1;
    int res_2;
};

typedef adjacency_list<vecS, vecS, directedS, SPPRC_Graph_Vert, SPPRC_Graph_Arc> SPPRCGraph;

struct Arc_edmonds {
    Arc_edmonds(int origin, int destine, double cost) : o(origin), d(destine), c(cost) {}
    int o, d;
    double c;
};

class Model {

public:
    typedef adjacency_list<vecS, vecS, directedS, property<vertex_index_t, int>, property<edge_weight_t, double>> BoostGraph;
    typedef graph_traits<BoostGraph>::edge_descriptor Edge;
    typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;

    typedef graph_traits<SPPRCGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<SPPRCGraph>::edge_descriptor edge_descriptor;

    // SPPRCGraph CshpGraph;
    vector<SPPRCGraph> cshpGraph;
    BoostGraph graphEdmonds;
    Graph *graph;

    property_map<BoostGraph, edge_weight_t>::type weightMap;
    property_map<BoostGraph, vertex_index_t>::type indexMap;
    vector<bool> z;
    vector<bool> noPath;
    vector<vector<bool>> y;
    vector<Arc_edmonds> arcsReturn;
    vector<vector<bool>> treeY;
    vector<vector<vector<bool>>> f;
    // vector<double> distance;
    // vector<size_t> parent;
    double objectiveValue = 0;

    void edmonds();
    
    void constrainedShortestpath(int k);

    void updateEdgePath(int i, int j, int k, double weight, bool increase);

    void updateEdgeBranching(int i,  int j, double weight);

    void insertPenalties(double penalty);

    Model(Graph *graph);

    void initialize();

    bool solve();

    bool isAcyclic();

    static bool compareArcs(Arc_edmonds a, Arc_edmonds b);
};

#endif //LAGRANGEANMST_MODEL_H
