
//
// Created by carlos on 28/05/19.
//

#include <iomanip>
#include "../headers/Model.h"

// data structures for shortest path problem with resource constraint
// ResourceContainer model
struct spp_spp_2_res_cont {
  spp_spp_2_res_cont(double c = 0, int r_1 = 0, int r_2 = 0) : cost(c), res_1(r_1), res_2(r_2) {}

  spp_spp_2_res_cont &operator=(const spp_spp_2_res_cont &other) {
    if (this == &other) return *this;
    this->~spp_spp_2_res_cont();
    new(this) spp_spp_2_res_cont(other);
    return *this;
  }

  double cost;
  int res_1;
  int res_2;
};

// data structures for shortest path problem with resource constraint
// ResourceContainer model
struct spp_spp_1_res_cont {
  spp_spp_1_res_cont(double c = 0, int r = 0) : 
    cost(c), res(r) {}

  spp_spp_1_res_cont &operator=(const spp_spp_1_res_cont &other) {
    if (this == &other) return *this;
    this->~spp_spp_1_res_cont();
    new(this) spp_spp_1_res_cont(other);
    return *this;
  }
  double cost;
  int res;
};

bool operator==(const spp_spp_2_res_cont &res_cont_1, const spp_spp_2_res_cont &res_cont_2) {
  return (res_cont_1.cost == res_cont_2.cost && res_cont_1.res_1 == res_cont_2.res_1 && res_cont_1.res_2 == res_cont_2.res_2);
}

bool operator==(const spp_spp_1_res_cont &res_cont_1, const spp_spp_1_res_cont &res_cont_2) {
  return (res_cont_1.cost == res_cont_2.cost && res_cont_1.res == res_cont_2.res);
}

bool operator<(const spp_spp_2_res_cont &res_cont_1, const spp_spp_2_res_cont &res_cont_2) {
  if (res_cont_1.cost > res_cont_2.cost) return false;
  if (res_cont_1.cost == res_cont_2.cost) return (res_cont_1.res_1 < res_cont_2.res_1 && res_cont_1.res_2 < res_cont_2.res_2);
  return true;
}

bool operator<(const spp_spp_1_res_cont &res_cont_1, const spp_spp_1_res_cont &res_cont_2) {
  if (res_cont_1.cost > res_cont_2.cost) return false;
  if (res_cont_1.cost == res_cont_2.cost) return (res_cont_1.res < res_cont_2.res);
  return true;
}

// ResourceExtensionFunction model
class ref_spprc_2_res {
public:
  inline bool operator()(const SPPRCGraph2Res &g, spp_spp_2_res_cont &new_cont, 
			 const spp_spp_2_res_cont &old_cont,
			 graph_traits<SPPRCGraph2Res>::edge_descriptor ed) const {

    const SPPRC_Graph_Arc_2_Res &arc_prop = get(edge_bundle, g)[ed];
    const SPPRC_Graph_Vert_2_Res &vert_prop = get(vertex_bundle, g)[target(ed, g)];
    new_cont.cost = old_cont.cost + arc_prop.cost;
    int &i_res_1 = new_cont.res_1;
    int &i_res_2 = new_cont.res_2;
    i_res_1 = old_cont.res_1 + arc_prop.res_1;
    i_res_2 = old_cont.res_2 + arc_prop.res_2;

    return (i_res_1 <= vert_prop.con_1 && i_res_2 <= vert_prop.con_2);
  }
};

// ResourceExtensionFunction model
class ref_spprc_1_res {
public:
  inline bool operator()(const SPPRCGraph1Res &g, spp_spp_1_res_cont &new_cont, 
			 const spp_spp_1_res_cont &old_cont,
			 graph_traits<SPPRCGraph1Res>::edge_descriptor ed) const {

    const SPPRC_Graph_Arc_1_Res &arc_prop = get(edge_bundle, g)[ed];
    const SPPRC_Graph_Vert_1_Res &vert_prop = get(vertex_bundle, g)[target(ed, g)];
    new_cont.cost = old_cont.cost + arc_prop.cost;
    int &i_res = new_cont.res;
    i_res = old_cont.res + arc_prop.res;

    return (i_res <= vert_prop.con);
  }
};

// DominanceFunction model
class dominance_spp_2_res {
public:
  inline bool operator()(const spp_spp_2_res_cont &res_cont_1, 
			 const spp_spp_2_res_cont &res_cont_2) const {
    return res_cont_1.cost <= res_cont_2.cost && 
      res_cont_1.res_1 <= res_cont_2.res_1 && 
      res_cont_1.res_2 <= res_cont_2.res_2;
  }
};
// end data structures for shortest path problem with time windows (spptw)

// DominanceFunction model
class dominance_spp_1_res {
public:
  inline bool operator()(const spp_spp_1_res_cont &res_cont_1, 
			 const spp_spp_1_res_cont &res_cont_2) const {
    return res_cont_1.cost <= res_cont_2.cost && res_cont_1.res <= res_cont_2.res;
  }
};

/*
  1 - Shortest Path + |V|-2 Arcs
  2 - Resource Constranined Shortest Path (one resource)
  2 - Resource Constranined Shortest Path (one resource)
  2 - Resource Constranined Shortest Path (two resources)
*/

Model::Model(Graph *graph, int relaxNum, bool heuristics) {
  Model::graph = graph;
  Model::heuristics = heuristics;
  Model::relaxNum = relaxNum;
  int n = graph->getN();

  if (heuristics) {
    heuristicGraph = BoostGraph();
    branchingEdges = vector<Arc*>();
    // Create vertex to heuristic graphs            
    graphEdmonds = BoostGraph(n);
    heuristicGraph = BoostGraph(n);
  }

  if (relaxNum == 1) createGraphRL1();
  else if (relaxNum == 2) createGraphRL2();
  else if (relaxNum == 3) createGraphRL3();
  else createGraphRL4();

  // Auxiliar vectors
  arcsReturn = vector<Arc_edmonds>();
  cout << "The graphs are created!" << endl;
}

void Model::createGraphRL1() {
  int n = graph->getN();
  shpGraph = vector<BoostGraph>(n);
  for (auto k : graph->DuS) shpGraph[k] = BoostGraph(n);

  for (int i = 0; i < n; i++) {
    for (auto arc : graph->arcs[i]) {
      // Heuristic Graphs
      if (heuristics) {
	add_edge(i, arc->getD(), 0, graphEdmonds);
	add_edge(i, arc->getD(), arc->getDelay(), heuristicGraph);
      }
      // Constrained SHP graph: Terminals
      for (auto k : graph->DuS)
	if (!graph->removedF[i][arc->getD()][k])
	  add_edge(i, arc->getD(), 0.0, shpGraph[k]);
    }
  }
}

void Model::createGraphRL2() {
  int n = graph->getN();
  cshpGraph1Res = vector<SPPRCGraph1Res>(n);
  shpGraph = vector<BoostGraph>(n);

  // Create vertex to the CSHP
  for (auto k : graph->terminals)
    for(int i = 0; i < n; i++)
      add_vertex(SPPRC_Graph_Vert_1_Res(i, graph->getParamDelay()), cshpGraph1Res[k]);
  // Create vertex to the SHP
  for (auto k : graph->DuS)
    shpGraph[k] = BoostGraph(n);        
  // Create the arcs
  int countEdges = 0;
  for (int i = 0; i < n; i++) {
    for (auto arc : graph->arcs[i]) {
      // Heuristic Graphs
      if (heuristics) {
	add_edge(i, arc->getD(), 0, graphEdmonds);
	add_edge(i, arc->getD(), arc->getDelay(), heuristicGraph);
      }
      // Constrained SHP graph: Terminals
      for (auto k : graph->terminals) {
	if (!graph->removedF[i][arc->getD()][k]) {
	  add_edge(i, arc->getD(), SPPRC_Graph_Arc_1_Res(countEdges++, 0, arc->getDelay()), cshpGraph1Res[k]);
	  add_edge(i, arc->getD(), 0.0, shpGraph[k]);
	}
      }
      for (auto q : graph->nonTerminals)
	add_edge(i, arc->getD(), 0.0, shpGraph[q]);
    }
  } 
}

void Model::createGraphRL3() {
  int n = graph->getN();
  cshpGraph1Res = vector<SPPRCGraph1Res>(n);
  shpGraph = vector<BoostGraph>(n);

  // Create vertex to the CSHP
  for (auto k : graph->terminals)
    for(int i = 0; i < n; i++)
      add_vertex(SPPRC_Graph_Vert_1_Res(i, graph->getParamJitter()), cshpGraph1Res[k]);
  // Create vertex to the SHP
  for (auto k : graph->DuS)
    shpGraph[k] = BoostGraph(n);        
  // Create the arcs
  int countEdges = 0;
  for (int i = 0; i < n; i++) {
    for (auto arc : graph->arcs[i]) {
      // Heuristic Graphs
      if (heuristics) {
	add_edge(i, arc->getD(), 0, graphEdmonds);
	add_edge(i, arc->getD(), arc->getDelay(), heuristicGraph);
      }
      // Constrained SHP graph: Terminals
      for (auto k : graph->terminals) {
	if (!graph->removedF[i][arc->getD()][k]) {
	  add_edge(i, arc->getD(), SPPRC_Graph_Arc_1_Res(countEdges++, 0, arc->getJitter()), cshpGraph1Res[k]);
	  add_edge(i, arc->getD(), 0.0, shpGraph[k]);
	}
      }
      for (auto q : graph->nonTerminals)
	add_edge(i, arc->getD(), 0.0, shpGraph[q]);
    }
  } 
}

void Model::createGraphRL4() {
  int n = graph->getN();
  cshpGraph2Res = vector<SPPRCGraph2Res>(n);
  shpGraph = vector<BoostGraph>(n);

  // Create vertex to the CSHP
  for (auto k : graph->terminals)
    for(int i = 0; i < n; i++)
      add_vertex(SPPRC_Graph_Vert_2_Res(i, graph->getParamDelay(), graph->getParamJitter()), cshpGraph2Res[k]);
  // Create vertex to the SHP
  for (auto k : graph->DuS)
    shpGraph[k] = BoostGraph(n);        
  // Create the arcs
  int countEdges = 0;
  for (int i = 0; i < n; i++) {
    for (auto arc : graph->arcs[i]) {
      // Heuristic Graphs
      if (heuristics) {
	add_edge(i, arc->getD(), 0, graphEdmonds);
	add_edge(i, arc->getD(), arc->getDelay(), heuristicGraph);
      }
      // Constrained SHP graph: Terminals
      for (auto k : graph->terminals) {
	if (!graph->removedF[i][arc->getD()][k]) {
	  add_edge(i, arc->getD(), SPPRC_Graph_Arc_2_Res(countEdges++, 0, arc->getDelay(), arc->getJitter()), cshpGraph2Res[k]);
	  add_edge(i, arc->getD(), 0.0, shpGraph[k]);
	}
      }
      for (auto q : graph->nonTerminals)
	add_edge(i, arc->getD(), 0.0, shpGraph[q]);
    }
  } 
}

void Model::initialize() {
  int n = graph->getN();
  y = vector<vector<bool>>(n, vector<bool>(n));
  treeY = vector<vector<bool>>(n, vector<bool>(n));
  z = vector<bool>(n);
  f = vector<vector<vector<bool>>>(n, vector<vector<bool>>(n, vector<bool>(n)));

  arcsReturn.erase(arcsReturn.begin(), arcsReturn.end());
  branchingEdges.erase(branchingEdges.begin(), branchingEdges.end());
}

bool Model::compareArcs(Arc_edmonds a, Arc_edmonds b) {
  return a.c < b.c;
}

double Model::arcsSelection() {
  double objective = 0;

  sort(arcsReturn.begin(), arcsReturn.end(), compareArcs);

  // Selecting |V| - 2 cause are one obrigatory arc
  for (int i = 0; i < graph->getN()-2; i++) {
    if (arcsReturn[i].o == graph->getRoot() && arcsReturn[i].d == 0) continue;
    objective += arcsReturn[i].c;
    y[arcsReturn[i].o][arcsReturn[i].d] = true;   
  }

  // Select the arcs (s, s')
  //  tie(ed, found) = edge(graph->getRoot(), 0, graphEdmonds);
  //cout << "e: " << endl;
  //getchar();
  //if (found) objective += boost::get(edge_weight_t(), graphEdmonds, ed);
  for (auto ac : arcsReturn)
    if (ac.o == graph->getRoot() && ac.d == 0) objective += ac.c;
  y[graph->getRoot()][0] = true;  
 
  return objective;
}

void Model::updateEdgeBranching(int i, int j, double weight) {
  if (heuristics) {
    Edge e;
    bool found;
    tie(e, found) = edge(i, j, graphEdmonds);
    if (found) {
      // Update the Edge to the edmonds graph and add to the set of arcs to be returned
      boost::put(edge_weight_t(), graphEdmonds, e, weight);
    }
  }
  //cout << i << " - " << j << " = " << weight << endl;
  arcsReturn.push_back(Arc_edmonds{i, j, weight});
}

void Model::edmonds() {
  int i, j;
  vector<Edge> branching;

  Vertex roots[2];
  roots[0] = graph->getRoot();
  roots[1] = graph->getRoot();

  edmonds_optimum_branching<false, true, true>(graphEdmonds,
					       indexMap,
					       weightMap,
					       roots,
					       roots+1,
					       back_inserter(branching));

  for (auto e : branching) {
    i = e.m_source, j = e.m_target;
    branchingEdges.push_back(new Arc(i, j, graph->getDelay(i, j), graph->getJitter(i, j), 0, 0));
    treeY[i][j] = true;
  }
}

double Model::shpTerminals(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  int actual, i, j, l, root = graph->getRoot();
  double minPath, minSP, objective;

  vector<double> distance(graph->getN(), (numeric_limits<int>::max)());
  vector<size_t> parent(graph->getN());

  for (i = 0; i < graph->getN(); ++i) parent[i] = i;
  distance[root] = 0;
    
  bool r = bellman_ford_shortest_paths(shpGraph[k], graph->getN(),
				       weight_map(weightMap).distance_map(&distance[0]).predecessor_map(&parent[0]));
    
  // If exists negative cycle, compute a cshp with values of bigM as resource
  //if (!r)
      //return makeModelRL1(k, multipliersDelay, multipliersJitter, multipliersVar, multipliersLeaf, multipliersRel);
  //else {
    actual = k;
    while(actual != root) {
      f[parent[actual]][actual][k] = true;
      actual = parent[actual];
    }
    return distance[k];
  //}
}

double Model::cshpTerminal1Res(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  int actual, indexMin, i, j, l, root = graph->getRoot();
  double minPath, minSP, objective;

  vector<double> distance(graph->getN(), (numeric_limits<int>::max)());
  vector<size_t> parent(graph->getN());

  for (i = 0; i < graph->getN(); ++i)
    parent[i] = i;
  distance[root] = 0;
    
  bool r = bellman_ford_shortest_paths(shpGraph[k], graph->getN(),
				       weight_map(weightMap).distance_map(&distance[0]).predecessor_map(&parent[0]));
    
  // If exists negative cycle, compute a cshp with values of bigM as resource
  if (!r) {
      cout << "Nagative" << endl;
   // if (relaxNum == 2) return makeModelRL2(k, multipliersDelay, multipliersVar, multipliersLeaf, multipliersRel);
   // else return makeModelRL3(k, multipliersJitter, multipliersVar, multipliersLeaf, multipliersRel);
  } else {
    vector<vector<graph_traits<SPPRCGraph1Res>::edge_descriptor>> opt_solutions;
    vector<spp_spp_1_res_cont> pareto_opt;
    
    SPPRC_Graph_Vert_1_Res &vert_prop = get(vertex_bundle, cshpGraph1Res[k])[k];
    if (relaxNum == 2) vert_prop.con = graph->getParamDelay();
    else vert_prop.con = graph->getParamJitter();

    r_c_shortest_paths(cshpGraph1Res[k],
		       get(&SPPRC_Graph_Vert_1_Res::num, cshpGraph1Res[k]),
		       get(&SPPRC_Graph_Arc_1_Res::num, cshpGraph1Res[k]),
		       root,
		       k,
		       opt_solutions,
		       pareto_opt,
		       spp_spp_1_res_cont(0, 0),
		       ref_spprc_1_res(),
		       dominance_spp_1_res(),
		       allocator<r_c_shortest_paths_label<SPPRCGraph1Res, spp_spp_1_res_cont >>(),
		       default_r_c_shortest_paths_visitor());

    if (!pareto_opt.empty()) {
      minPath = pareto_opt[0].cost;
      indexMin = 0;
      for (i = 0; i < opt_solutions.size(); i++) {
	if (pareto_opt[i].cost < minPath) {
	  minPath = pareto_opt[i].cost;
	  indexMin = i;
	}
      }

      if (minPath < distance[k]) {
	for (j = static_cast<int>(opt_solutions[indexMin].size())-1; j >= 0; --j) {
	  i = source(opt_solutions[indexMin][j], cshpGraph1Res[k]);
	  l = target(opt_solutions[indexMin][j], cshpGraph1Res[k]);
	  f[i][l][k] = true;
	}
	
	return minPath;
      } else {
	actual = k;
	while(actual != root) {
	  f[parent[actual]][actual][k] = true;
	  actual = parent[actual];
	}
	return distance[k];
      }
    } else {
      z[k] = true;
      actual = k;
      while(actual != root) {
	f[parent[actual]][actual][k] = true;
	actual = parent[actual];
      }
      return distance[k];
    }
  }
}

double Model::cshpTerminal2Res(int k, vector<vector<double>> &multipliersVar, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  int actual, indexMin, i, j, l, root = graph->getRoot();
  double minPath, minSP, objective;
  vector<pair<int, int>> firstPath = vector<pair<int, int>>();

  vector<double> distance(graph->getN(), (numeric_limits<int>::max)());
  vector<size_t> parent(graph->getN());

  for (i = 0; i < graph->getN(); ++i)
    parent[i] = i;
  distance[root] = 0;
    
  bool r = bellman_ford_shortest_paths(shpGraph[k], graph->getN(),
				       weight_map(weightMap).distance_map(&distance[0]).predecessor_map(&parent[0]));
  
  // If exists negative cycle, compute a cshp with values of bigM as resource
  if (!r) {
    //return makeModelRL4(k, multipliersVar, multipliersLeaf, multipliersRel);
  } else {
    vector<vector<graph_traits<SPPRCGraph2Res>::edge_descriptor>> opt_solutions;
    vector<spp_spp_2_res_cont> pareto_opt;
    
    SPPRC_Graph_Vert_2_Res &vert_prop = get(vertex_bundle, cshpGraph2Res[k])[k];
    vert_prop.con_1 = graph->getParamDelay();
    vert_prop.con_2 = graph->getParamJitter();

    r_c_shortest_paths(cshpGraph2Res[k],
		       get(&SPPRC_Graph_Vert_2_Res::num, cshpGraph2Res[k]),
		       get(&SPPRC_Graph_Arc_2_Res::num, cshpGraph2Res[k]),
		       root,
		       k,
		       opt_solutions,
		       pareto_opt,
		       spp_spp_2_res_cont(0, 0, 0),
		       ref_spprc_2_res(),
		       dominance_spp_2_res(),
		       allocator<r_c_shortest_paths_label<SPPRCGraph2Res, spp_spp_2_res_cont >>(),
		       default_r_c_shortest_paths_visitor());

    if (!pareto_opt.empty()) {
      minPath = pareto_opt[0].cost;
      indexMin = 0;
      for (i = 0; i < opt_solutions.size(); i++) {
	if (pareto_opt[i].cost < minPath) {
	  minPath = pareto_opt[i].cost;
	  indexMin = i;
	}
      }

      if (minPath < distance[k]+1) {
	for (j = static_cast<int>(opt_solutions[indexMin].size())-1; j >= 0; --j) {
	  i = source(opt_solutions[indexMin][j], cshpGraph2Res[k]);
	  l = target(opt_solutions[indexMin][j], cshpGraph2Res[k]);
	  f[i][l][k] = true;
	}
	return minPath;
      } else {
	actual = k;
    z[k] = true;
	while(actual != root) {
	  f[parent[actual]][actual][k] = true;
	  actual = parent[actual];
	}
	return distance[k];
      }
    } else {
      z[k] = true;
      actual = k;
      while(actual != root) {
	f[parent[actual]][actual][k] = true;
	actual = parent[actual];
      }
      return distance[k];
    }
  }
}

double Model::shpNonTerminal(int q) {
  int n = graph->getN(), root = graph->getRoot(), actual;
  property_map<BoostGraph, edge_weight_t>::type weightMap = get(edge_weight, shpGraph[q]);
  vector<Vertex> predecessors = vector<Vertex>(n);
  vector<double> distance = vector<double>(n);

  dijkstra_shortest_paths(shpGraph[q], root, predecessor_map(
							     make_iterator_property_map(predecessors.begin(), get(vertex_index, shpGraph[q]))).distance_map(
																			    make_iterator_property_map(distance.begin(), get(vertex_index, shpGraph[q]))));

  actual = q;
  while(actual != root) {
    f[predecessors[actual]][actual][q] = true;
    actual = predecessors[actual];
  }

  return distance[q];
}

void Model::updateEdgePath(int i, int j, int k, bool isTerminal, double weight) {
  edge_descriptor_1_res ed1Res;
  edge_descriptor_2_res ed2Res;
  Edge ed;
  bool found;

  if (isTerminal) {
    if (relaxNum <= 3) tie(ed1Res, found) = edge(i,  j, cshpGraph1Res[k]);
    else tie(ed2Res, found) = edge(i,  j, cshpGraph2Res[k]);

    if (found) {
      if (relaxNum <= 3) {
	SPPRC_Graph_Arc_1_Res &arc_prop = get(edge_bundle, cshpGraph1Res[k])[ed1Res];
	arc_prop.cost = weight;
      } else {
	SPPRC_Graph_Arc_2_Res &arc_prop = get(edge_bundle, cshpGraph2Res[k])[ed2Res];
	arc_prop.cost = weight;
      }
    }
  }

  tie(ed, found) = boost::edge(i, j, shpGraph[k]);
  if (found) boost::put(edge_weight_t(), shpGraph[k], ed, weight);
}

double Model::getEdgeShp(int i, int j, int k) {
  Edge ed;
  bool found;
  double result = 0;
    
  tie(ed, found) = edge(i, j, shpGraph[k]);
  if (found) result = boost::get(edge_weight_t(), shpGraph[k], ed);
    
  return result;
}

double Model::getEdgeCshp(int i, int j, int k) {
  edge_descriptor_1_res ed1Res;
  edge_descriptor_2_res ed2Res;
  bool found;
  double cost = 0;

  if (relaxNum <= 3) tie(ed1Res, found) = edge(i,  j, cshpGraph1Res[k]);
  else tie(ed2Res, found) = edge(i,  j, cshpGraph2Res[k]);

  if (found) {
    if (relaxNum <= 3) {
      SPPRC_Graph_Arc_1_Res &arc_prop = get(edge_bundle, cshpGraph1Res[k])[ed1Res];
      cost = arc_prop.cost;
    } else {
      SPPRC_Graph_Arc_2_Res &arc_prop = get(edge_bundle, cshpGraph2Res[k])[ed2Res];
      cost = arc_prop.cost;
    }
  }

  return cost;
}

bool Model::isAcyclic() {
  int n = graph->getN();
  BoostGraph G(n);

  for (int i = 0; i < graph->getN()-1; i++)
    add_edge(arcsReturn[i].o, arcsReturn[i].d, 1, G);

  vector<int> component(n), discover_time(n);
  vector<default_color_type> color(n);
  vector<Vertex> rootG(n);

  int num = strong_components(G, make_iterator_property_map(component.begin(), get(vertex_index, G)),
			      root_map(make_iterator_property_map(rootG.begin(), get(vertex_index, G))).
			      color_map(make_iterator_property_map(color.begin(), get(vertex_index, G))).
			      discover_time_map(make_iterator_property_map(discover_time.begin(), get(vertex_index, G))));

  return !(num != n);
}

int Model::initialHeuristic() {
  int n = graph->getN();
  property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, heuristicGraph);
  vector<Vertex> predecessors = vector<Vertex>(n);
  vector<int> distance = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n);
  vector<int> jitterDistance = vector<int>(n);

  dijkstra_shortest_paths(heuristicGraph, graph->getRoot(), predecessor_map(
									    make_iterator_property_map(predecessors.begin(), get(vertex_index, heuristicGraph))).distance_map(
																					      make_iterator_property_map(distance.begin(), get(vertex_index, heuristicGraph))));

  double meanValue = 0;
  int countTerm = 0, actual, jitter, count = 0;

  for (auto k : graph->terminals) 
    if (distance[k] > graph->getParamDelay()) {
      countTerm++;
      notAttended[k] = true;
    } else meanValue += distance[k];

  for (auto t : graph->terminals) {
    actual = t, jitter = 0;
    while (actual != graph->getRoot()) {
      jitter += graph->getJitter(predecessors[actual], actual);
      actual = predecessors[actual];
    }
    if (jitter > graph->getParamJitter()) 
      notAttended[t] = true;
    jitterDistance[t] = jitter;
  }

  meanValue = meanValue / double(graph->terminals.size() - countTerm);
  int lessThanAvg = 0, greaThanAvg = 0;
    
  for (auto k : graph->terminals) {
    if (!notAttended[k]) {
      if (distance[k] <= meanValue) lessThanAvg++;
      else greaThanAvg++;
    }
  }

  for (auto k : graph->terminals) {
    if (!notAttended[k]) {
      for (auto l : graph->terminals) {
	if (l != k && !notAttended[l]) {
	  if (distance[k] - distance[l] > graph->getParamVariation()) {
	    if (lessThanAvg < greaThanAvg) {
	      if (distance[k] < distance[l]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    } else {
	      if (distance[k] > distance[l]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    }
	  }
	}
      }
    }
  }

  for (auto t : graph->terminals) 
    if (notAttended[t]) count++;  
  return count;
}

int Model::subgradientHeuristic() {
  cout << "heuristic " << endl;
  int i, j, n = graph->getN(), bestCand, root = graph->getRoot();
  vector<int> delayPaths = vector<int>(n), jitterPaths = vector<int>(n),
    predecessors = vector<int>(n), delayPathAux = vector<int>(n),
    jitterPathAux = vector<int>(n), predecessorsAux = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n), notAttendedAux = vector<bool>(n);
  vector<int> changed = vector<int>();

  predecessors[root] = root;
  for (auto arc : branchingEdges) {
    predecessors[arc->getD()] = arc->getO();
  }

  random_shuffle(graph->terminals.begin(), graph->terminals.end());

  int actual, jitter, delay, count = 0;
  for (auto k : graph->DuS) {
    actual = k, jitter = 0, delay = 0;
    while (actual != root) {
      jitter += graph->getJitter(predecessors[actual], actual),
	delay += graph->getDelay(predecessors[actual], actual);
      actual = predecessors[actual];
    }
    delayPaths[k] = delay, jitterPaths[k] = jitter;
    if (delay > graph->getParamDelay() || jitter > graph->getParamJitter())
        notAttended[k] = true;
  }
 
  // Select a terminal to fix as attend
  int diffDelay, diffJitter;
  bool canMove;
  int selected = -1, losts;

  for (auto k : graph->terminals) {
    if (!notAttended[k]) {
      selected = k;
      break;
    }
  }

  if (selected == -1) 
    return graph->terminals.size();

  // get the values of this path
  for (auto k : graph->terminals) {
    if (k != selected) {
      bestCand = -1, delay = graph->getParamDelay()+1, jitter = graph->getParamJitter()+1; 
            
      if (delayPaths[k] < (delayPaths[selected] - graph->getParamVariation()) ||
	  delayPaths[k] > (delayPaths[selected] + graph->getParamVariation())) {
                
	for (auto arc : graph->arcs[k]) {
	  i = arc->getD();
	  // Get the best candidate to move
	  if (i != predecessors[k] && k != predecessors[i]) {
	    if (delayPaths[i] + arc->getDelay() >= (delayPaths[selected] - graph->getParamVariation()) &&
		delayPaths[i] + arc->getDelay() <= (delayPaths[selected] + graph->getParamVariation()) &&
		delayPaths[i] + arc->getDelay() <= graph->getParamDelay() && 
		jitterPaths[i] + arc->getJitter() <= graph->getParamJitter()) {

	      bestCand = i;
	      delay = delayPaths[i] + arc->getDelay();
	      jitter = jitterPaths[i] + arc->getJitter();
	      break;
	    }
	  }
	}

	if (bestCand != -1) {
	  // create the temporary vectors
	  canMove = true;
	  losts = 0;
	  notAttended[k] = false;
	  vector<vector<int>> sub = vector<vector<int>>(n, vector<int>());
	  for (i = 0; i < n; i++) {
	    delayPathAux[i] = delayPaths[i], jitterPathAux[i] = jitterPaths[i], 
	      predecessorsAux[i] = predecessors[i], notAttendedAux[i] = notAttended[i];
	    sub[predecessors[i]].push_back(i);
	  }

	  // Evaluate the move                    
	  diffDelay = delay - delayPathAux[k], diffJitter = jitter - jitterPathAux[k];
	  predecessorsAux[k] = bestCand;
	  delayPathAux[k] = delay, jitterPathAux[k] = jitter;

	  changed.erase(changed.begin(), changed.end());
	  changed.push_back(k);
	  
	  while (!changed.empty()) {
	    actual = changed.back();
	    changed.pop_back();

	    for (int j : sub[actual]) {
	      delayPathAux[j] += diffDelay, jitterPathAux[j] += diffJitter;
	      if (delayPathAux[j] > graph->getParamDelay() || jitterPathAux[j] > graph->getParamJitter()) {
		if (!notAttendedAux[j]) {
		  notAttendedAux[j] = true;
		  losts++;
		}
	      } else notAttendedAux[j] = false;
	      if (losts >= 2 || (j == selected && notAttendedAux[j])) {
		  canMove = false;
		  changed.erase(changed.begin(), changed.end());
		  break;
	      } else changed.push_back(j);	      
	    }
	  }
	  if (canMove) {
	    for (i = 0; i < n; i++) {
	      delayPaths[i] = delayPathAux[i], jitterPaths[i] = jitterPathAux[i], 
		predecessors[i] = predecessorsAux[i], notAttended[i] = notAttendedAux[i];
	    }
	  }
	}  
      } 
    }
  }

  vector<int> delays = vector<int>();

  for (auto k : graph->terminals)
    if (delayPaths[k] <= graph->getParamDelay() && jitterPaths[k] <= graph->getParamJitter())
      delays.push_back(delayPaths[k]);
  
  sort(delays.begin(), delays.end());

  int actMax = 0, p2 = 0, maxi = 0, fd;
  for (i = 0; i < delays.size(); i++) {
    fd = delays[i];

    while (delays[p2] <= fd + graph->getParamVariation() && p2 < delays.size()) p2++;

    if ((p2 - i) > maxi) maxi = (p2 - i);
  }
  return graph->terminals.size() - maxi;
}

double Model::penalty1Rl (vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar) {
  double penalty = 0, coefZ, objective = 0;
  int bigMK, bigML;

  for (auto k : graph->terminals) {
    coefZ = -(multipliersDelay[k] * graph->getParamDelay()) - (multipliersJitter[k] * graph->getParamJitter());
    penalty += coefZ;
     
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                
	coefZ -= (multipliersVar[k][l] * bigMK + multipliersVar[l][k] * bigML);
	penalty -= (multipliersVar[k][l] * graph->getParamVariation());
      }
    }

    if (1 + coefZ < 0 || z[k] || graph->noPath[k]) {
      z[k] = true;
      objective += 1 + coefZ;
    }
  }

  return (objective + penalty);
}

double Model::penalty2Rl (vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar) {
  double penalty = 0, coefZ, objective = 0;
  int bigMK, bigML;

  for (auto k : graph->terminals) {
    coefZ = -(multipliersJitter[k] * graph->getParamJitter());
    penalty += coefZ;
     
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                
	coefZ -= (multipliersVar[k][l] * bigMK + multipliersVar[l][k] * bigML);
	penalty -= (multipliersVar[k][l] * graph->getParamVariation());
      }
    }

    if (1 + coefZ < 0 || z[k] || graph->noPath[k]) {
      z[k] = true;
      objective += 1 + coefZ;
    }
  }
    
  return (objective + penalty);
}

double Model::penalty3Rl (vector<double> &multipliersDelay, vector<vector<double>> &multipliersVar) {
  double penalty = 0, coefZ, objective = 0;
  int bigMK, bigML;

  for (auto k : graph->terminals) {
    coefZ = -(multipliersDelay[k] * graph->getParamDelay());
    penalty += coefZ;
     
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                
	coefZ -= (multipliersVar[k][l] * bigMK + multipliersVar[l][k] * bigML);
	penalty -= (multipliersVar[k][l] * graph->getParamVariation());
      }
    }

    if (1 + coefZ < 0|| z[k] || graph->noPath[k]) {
      z[k] = true;
      objective += 1 + coefZ;
    }
  }
    
  return (objective + penalty);
}

double Model::penalty4Rl (vector<vector<double>> &multipliersVar) {
  double penalty = 0, coefZ, objective = 0;
  int bigMK, bigML;

  for (auto k : graph->terminals) {
    coefZ = 0;     
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                
	coefZ -= (multipliersVar[k][l] * bigMK + multipliersVar[l][k] * bigML);
	penalty -= (multipliersVar[k][l] * graph->getParamVariation());
      }
    }

    if (1 + coefZ < 0 || z[k] || graph->noPath[k]) {
      z[k] = true;
      objective += 1 + coefZ;
    }
  }
  return (objective + penalty);
}

void Model::updateArbCosts(vector<vector<vector<double>>> &multipliersRel) { 
  int i, j, o, d, n = graph->getN();
  double sumMultipliers;

  for (i = 0; i < n; i++) {
    for (auto *arc : graph->arcs[i]) {
      j = arc->getD();
      sumMultipliers = 0;
      for (auto k : graph->DuS)
	    sumMultipliers -= multipliersRel[i][j][k];
      //      cout << i << " - " << j << " = " << sumMultipliers << endl;
      updateEdgeBranching(i, j, sumMultipliers);
    }
  }
  cout << "Arb updated" << endl;
}

void Model::updatePathCostsNT(int q, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  double edgeCost = 0;

  for (auto *arc : graph->arcs[0])
    updateEdgePath(0, arc->getD(), q, false, multipliersRel[0][arc->getD()][q] + multipliersLeaf[arc->getD()][q]);

  for (int i = 0; i < graph->getN(); i++)
    for (auto *arc : graph->arcs[i]) 
      updateEdgePath(i, arc->getD(), q, false, multipliersRel[i][arc->getD()][q]);
}

void Model::updatePathCostsTerminals(int k, vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  double edgeCostAux = 0;
  int i, j, o, d, n = graph->getN();

  // Default for all relaxations
  for (auto l : graph->terminals)
    if (l != k) edgeCostAux += (multipliersVar[k][l] - multipliersVar[l][k]);

  // SHP + Selection
  if (relaxNum == 1) {
    for (auto *arc : graph->arcs[0]) {
      j = arc->getD();
      if (j != k)
	updateEdgePath(0, j, k, true, (edgeCostAux + multipliersDelay[k]) * arc->getDelay() +
		       multipliersJitter[k] * arc->getJitter() + multipliersRel[0][j][k] + multipliersLeaf[j][k]);
    }

    for (i = 1; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	updateEdgePath(i, j, k, true, (edgeCostAux + multipliersDelay[k]) * arc->getDelay() + 
		       multipliersJitter[k] * arc->getJitter() + multipliersRel[i][j][k]);
      }
    }
  } else if (relaxNum == 2) {
    for (auto *arc : graph->arcs[0]) {
      j = arc->getD();
      if (j != k)
	updateEdgePath(0, j, k, true, edgeCostAux * arc->getDelay() + 
		       multipliersJitter[k] * arc->getJitter() + multipliersRel[0][j][k] + multipliersLeaf[j][k]);
    }

    for (i = 1; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	updateEdgePath(i, j, k, true, edgeCostAux * arc->getDelay() + multipliersJitter[k] * arc->getJitter() + multipliersRel[i][j][k]);
      }
    }
  } else if (relaxNum == 3) {
    for (auto *arc : graph->arcs[0]) {
      j = arc->getD();
      if (j != k)
	updateEdgePath(0, j, k, true, (edgeCostAux + multipliersDelay[k]) * arc->getDelay() + 
		       multipliersRel[0][j][k] + multipliersLeaf[j][k]);
    }

    for (i = 1; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	updateEdgePath(i, j, k, true, (edgeCostAux + multipliersDelay[k]) * arc->getDelay() + multipliersRel[i][j][k]);
      }
    }
  } else {
    for (auto *arc : graph->arcs[0]) {
      j = arc->getD();
       if (j != k) updateEdgePath(0, j, k, true, multipliersRel[0][j][k] + multipliersLeaf[j][k]);
    }

    for (i = 1; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
        j = arc->getD();
        updateEdgePath(i, j, k, true, multipliersRel[i][j][k]);
        //cout << i << " - " << j << " - " << k << " = " << multipliersRel[i][j][k] << endl;
      }
    }
  }
  //getchar();
}

bool Model::solve(vector<double> &multipliersDelay, vector<double> &multipliersJitter, vector<vector<double>> &multipliersVar, vector<vector<double>> &multipliersLeaf, vector<vector<vector<double>>> &multipliersRel) {
  initialize();
  objectiveFunction = 0;
  // Terminals 
  for (auto k : graph->terminals) {
    if (graph->noPath[k]) {
      f[graph->getRoot()][0][k] = f[0][k][k] = true;
      objectiveFunction += (getEdgeShp(graph->getRoot(), 0, k) + getEdgeShp(0, k, k));
      z[k] = true;
    } else {
      // Path to terminals
      updatePathCostsTerminals(k, multipliersDelay, multipliersJitter, multipliersVar, multipliersLeaf, multipliersRel);
      if (relaxNum == 1) objectiveFunction += shpTerminals(k, multipliersDelay, multipliersJitter, multipliersVar, multipliersLeaf, multipliersRel);
      else if (relaxNum <= 3) objectiveFunction += cshpTerminal1Res(k, multipliersDelay, multipliersJitter, multipliersVar, multipliersLeaf, multipliersRel);
      else objectiveFunction += cshpTerminal2Res(k, multipliersVar, multipliersLeaf, multipliersRel);
    }
  }
  //cout << "CSHP: " << objectiveFunction << endl;
  // Non-terminals
  for (auto k : graph->nonTerminals) {
    if (graph->removed[k]) {
      f[graph->getRoot()][0][k] = f[0][k][k] = true;
      objectiveFunction += (getEdgeShp(graph->getRoot(), 0, k) + getEdgeShp(0, k, k));
    } else {
      // Path to non-terminals
      updatePathCostsNT(k, multipliersLeaf, multipliersRel);
      objectiveFunction += shpNonTerminal(k);
    }
  }
  //cout << "SHP: " << objectiveFunction << endl;

  // Tree and Arc selection
  updateArbCosts(multipliersRel);
  objectiveFunction += arcsSelection();
  
  // Heuristics
  if (heuristics) {
    edmonds();
    heuristicObj = subgradientHeuristic();
  }
  
  //cout << "Arb: " << objectiveFunction << endl;
  // Penalties
  if (relaxNum == 1) objectiveFunction += penalty1Rl(multipliersDelay, multipliersJitter, multipliersVar);
  else if (relaxNum == 2) objectiveFunction += penalty2Rl(multipliersJitter, multipliersVar);
  else if (relaxNum == 3) objectiveFunction += penalty3Rl(multipliersDelay, multipliersVar);
  else {
    for (auto k : graph->terminals)
      if (z[k]) {
          objectiveFunction += 1;
      }
  }
  //cout << "PPL: " << objectiveFunction << endl;
  return true;
}

int Model::getOriginalObj() {
  int objective = 0;
  for (auto k : graph->terminals)
    if (z[k]) objective++;
  return objective;
}

double Model::getObj() {
  return objectiveFunction;
}

int Model::getHeuristicObj() {
  return heuristicObj;
}
