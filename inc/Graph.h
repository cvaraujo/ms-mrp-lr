//
// Created by carlos on 27/05/19.
//

#ifndef LAGRANGEANMST_GRAPH_H
#define LAGRANGEANMST_GRAPH_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <bits/ios_base.h>
#include <algorithm>
#include <fstream>
#include "Arc.h"
#include "edmonds_optimum_branching.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/config.hpp"
#include "boost/graph/r_c_shortest_paths.hpp"

#ifdef BOOST_MSVC
#  pragma warning(disable: 4267)
#endif

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int>> BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor EdgeDescriptor;
typedef graph_traits<BoostGraph>::vertex_descriptor VertexDescriptor;

class Graph {
    int n, m, root, cntRemoved, paramDelay, paramJitter, paramVariation, paramBandwidth, bigMDelay = 0, bigMJitter = 0;

    typedef graph_traits<BoostGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<BoostGraph>::edge_descriptor edge_descriptor;

    vector<vertex_descriptor> predecessors;
    vector<int> distance;

public:
    BoostGraph preProcessing;
    vector<vector<Arc *>> arcs;
    vector<int> terminals, nonTerminals, DuS, delayVector, jitterVector;
    vector<bool> removed;
    vector<vector<vector<bool>>> removedF;

    Graph(string instance, string param, string outputName);

    void SAE();

    void graphReduction();

    void showGraph();

    int getN() const;

    void setN(int n);

    int getM() const;

    void setM(int m);

    int getParamDelay() const;

    void setParamDelay(int paramDelay);

    int getParamJitter() const;

    void setParamJitter(int paramJitter);

    int getParamVariation() const;

    void setParamVariation(int paramVariation);

    int getParamBandwidth() const;

    void setParamBandwidth(int paramBandwidth);

    int getRoot() const;

    void setRoot(int root);

    int getBigMDelay();

    int getBigMJitter();

    int getShpTerminal(int k);

    int getNAfterRemoved();
};


#endif //LAGRANGEANMST_GRAPH_H
