//
// Created by carlos on 28/05/19.
//

#ifndef LAGRANGEANMST_MODEL_H
#define LAGRANGEANMST_MODEL_H

#include "Graph.h"
#include "gurobi_c++.h"

using namespace std;
using namespace boost;

struct SPPRC_Graph_Vert_2_Res {
  SPPRC_Graph_Vert_2_Res(int n = 0, int c_1 = 0, int c_2 = 0) : 
    num(n), con_1(c_1), con_2(c_2) {}

  int num;
  // Resource consumed
  int con_1;
  int con_2;
};

struct SPPRC_Graph_Arc_2_Res {
  SPPRC_Graph_Arc_2_Res(int n = 0, double c = 0, int r_1 = 0, int r_2 = 0) : 
    num(n), cost(c), res_1(r_1), res_2(r_2) {}

  int num;
  // traversal cost
  double cost;
  // traversal resource
  int res_1;
  int res_2;
};

struct SPPRC_Graph_Vert_1_Res {
  SPPRC_Graph_Vert_1_Res(int n = 0, int c = 0) : 
    num(n), con(c) {}

  int num;
  // Resource consumed
  int con;
};

struct SPPRC_Graph_Arc_1_Res {
  SPPRC_Graph_Arc_1_Res(int n = 0, double c = 0, int r = 0) : 
    num(n), cost(c), res(r) {}

  int num;
  // traversal cost
  double cost;
  // traversal resource
  int res;
};

typedef adjacency_list<vecS, vecS, directedS, SPPRC_Graph_Vert_1_Res, SPPRC_Graph_Arc_1_Res> SPPRCGraph1Res;
typedef adjacency_list<vecS, vecS, directedS, SPPRC_Graph_Vert_2_Res, SPPRC_Graph_Arc_2_Res> SPPRCGraph2Res;

struct Arc_edmonds {
  Arc_edmonds(int origin, int destine, double cost) : 
    o(origin), d(destine), c(cost) {}
  int o, d;
  double c;
};

class Model {
  bool heuristics;
  int relaxNum, originalObj, heuristicObj = -1;
  double objectiveFunction = 0;

  void createGraphRL1();
  void createGraphRL2();
  void createGraphRL3();
  void createGraphRL4();
  double shpTerminals(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);
  double cshpTerminal1Res(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);
  double cshpTerminal2Res(int k, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);
  double shpNonTerminal(int q);
  double getEdgeShp(int i, int j, int k);
  double getEdgeCshp(int i, int j, int k);

  double penalty1Rl (vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar);
  double penalty2Rl (vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar);
  double penalty3Rl (vector<double> &multipliersDelay, vector<vector<double>> &multipliersVar);
  double penalty4Rl (vector<vector<double>> &multipliersVar);

  void updateArbCosts(vector<vector<vector<double>>> &multipliersRel);
  void updatePathCostsNT(int q,  vector<vector<vector<double>>> &multipliersRel);
  void updatePathCostsTerminals(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, 
				vector<vector<double>> &multipliersVar,  
				vector<vector<vector<double>>> &multipliersRel);
  int subgradientHeuristic();

public:
  // Graph for simple shortest path
  typedef adjacency_list<vecS, vecS, directedS, property<vertex_index_t, int>, property<edge_weight_t, double>> BoostGraph;

  // Structs to get the vertex and edges
  typedef graph_traits<BoostGraph>::edge_descriptor Edge;
  typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;
  typedef graph_traits<SPPRCGraph1Res>::vertex_descriptor vertex_descriptor_1_res;
  typedef graph_traits<SPPRCGraph1Res>::edge_descriptor edge_descriptor_1_res;
  typedef graph_traits<SPPRCGraph2Res>::vertex_descriptor vertex_descriptor_2_res;
  typedef graph_traits<SPPRCGraph2Res>::edge_descriptor edge_descriptor_2_res;

  // Shortest Path with resource constraint 1 and 2 resources
  vector<SPPRCGraph1Res> cshpGraph1Res;
  vector<SPPRCGraph2Res> cshpGraph2Res;
    
  // Shorteste path graph
  vector<BoostGraph> shpGraph;

  // Heuristics graphs
  BoostGraph graphEdmonds;
  BoostGraph heuristicGraph;

  // Original graph (Backup pointer)
  Graph *graph;

  // Auxiliar maps from Boost Library
  property_map<BoostGraph, edge_weight_t>::type weightMap;
  property_map<BoostGraph, vertex_index_t>::type indexMap;

  // Variables from model a
  vector<bool> z;
  vector<vector<bool>> y, treeY;
  vector<vector<vector<bool>>> f;

  // auxiliar vector 
  vector<Arc_edmonds> arcsReturn;
  vector<Arc*> branchingEdges;
    
  Model(Graph *graph, int relaxNum, bool heuristics);

  double arcsSelection();

  void edmonds();
    
  void constrainedShortestpath(int k, bool isNonTerminal);

  void updateEdgePath(int i, int j, int k, bool isTerminal, double weight);

  double getEdgePath(int i, int j, int k);

  void updateEdgeBranching(int i,  int j, double weight);

  void insertPenalties(double penalty);

  void initialize();

  bool solve(vector<double> &multipliersDelay, vector<double> &multipliersJitter, 
	     vector<vector<double>> &multipliersVar,  
	     vector<vector<vector<double>>> &multipliersRel);

  double makeModelRL1(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);
  
  double makeModelRL2(int k, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);

  double makeModelRL3(int k, vector<double> &multipliersDelay, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);

  double makeModelRL4(int k, vector<vector<double>> &multipliersVar,  vector<vector<vector<double>>> &multipliersRel);
  
  int getOriginalObj();
    
  double getObj();

  int getHeuristicObj();
    
  int initialHeuristic();

  bool isAcyclic();

  static bool compareArcs(Arc_edmonds a, Arc_edmonds b);
};

#endif //LAGRANGEANMST_MODEL_H
