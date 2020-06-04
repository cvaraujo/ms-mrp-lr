//
// Created by carlos on 28/05/19.
//

#include <iomanip>
#include "Model.h"

// data structures for shortest path problem with resource constraint
// ResourceContainer model
struct spp_spp_res_cont {
    spp_spp_res_cont(double c = 0, int r_1 = 0, int r_2 = 0) : cost(c), res_1(r_1), res_2(r_2) {}

    spp_spp_res_cont &operator=(const spp_spp_res_cont &other) {
        if (this == &other) return *this;
        this->~spp_spp_res_cont();
        new(this) spp_spp_res_cont(other);
        return *this;
    }

    double cost;
    int res_1;
    int res_2;
};

bool operator==(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    return (res_cont_1.cost == res_cont_2.cost && res_cont_1.res_1 == res_cont_2.res_1 && res_cont_1.res_2 == res_cont_2.res_2);
}

bool operator<(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    if (res_cont_1.cost > res_cont_2.cost) return false;
    if (res_cont_1.cost == res_cont_2.cost) return (res_cont_1.res_1 < res_cont_2.res_1 && res_cont_1.res_2 < res_cont_2.res_2);
    return true;
}

// ResourceExtensionFunction model
class ref_spprc {
public:
    inline bool operator()(const SPPRCGraph &g, spp_spp_res_cont &new_cont, const spp_spp_res_cont &old_cont,
                           graph_traits<SPPRCGraph>::edge_descriptor ed) const {

        const SPPRC_Graph_Arc &arc_prop = get(edge_bundle, g)[ed];
        const SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.cost = old_cont.cost + arc_prop.cost;
        int &i_res_1 = new_cont.res_1;
        int &i_res_2 = new_cont.res_2;
        i_res_1 = old_cont.res_1 + arc_prop.res_1;
        i_res_2 = old_cont.res_2 + arc_prop.res_2;

        return (i_res_1 <= vert_prop.con_1 && i_res_2 <= vert_prop.con_2);
    }
};

// DominanceFunction model
class dominance_spptw {
public:
    inline bool operator()(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) const {
        return res_cont_1.cost <= res_cont_2.cost && res_cont_1.res_1 <= res_cont_2.res_1 && res_cont_1.res_2 <= res_cont_2.res_2;
    }
};
// end data structures for shortest path problem with time windows (spptw)

Model::Model(Graph *graph) {
    this->graph = graph;
    int n = graph->getN();
    cshpGraph = vector<SPPRCGraph>(n);
    vector<BoostGraph> aux = vector<BoostGraph>(n);
    graphEdmonds = BoostGraph();
    heuristicGraph = BoostGraph();
    noPath = vector<bool>(n);
    used = vector<bool>(n);
    branchingEdges = vector<Arc*>();

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            add_vertex(j, aux[i]);

    for (int i = 0; i < n; i++)
        for (auto arc : graph->arcs[i]) 
            for (auto k : graph->terminals)
                if (!graph->removedF[i][arc->getD()][k])
                    add_edge(i, arc->getD(), 1, aux[k]);
    
    for (auto k : graph->terminals) {
        noPath[k] = false;
        property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, aux[k]);
        vector<Vertex> predecessors = vector<Vertex>(n);
        vector<int> distance = vector<int>(n);

        dijkstra_shortest_paths(aux[k], graph->getRoot(), predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, aux[k]))).distance_map(
            make_iterator_property_map(distance.begin(), get(vertex_index, aux[k]))));

        if (distance[k] >= numeric_limits<int>::max()) {
            noPath[k] = true;
            cout << k+1 << endl;
        }
    }

    for (auto k : graph->DuS)
        if (!noPath[k] && !graph->removed[k]) 
            for(int i = 0; i < n; i++)
                add_vertex(SPPRC_Graph_Vert(i, graph->getParamDelay(), graph->getParamJitter()), cshpGraph[k]);
                // add_vertex(SPPRC_Graph_Vert(i, graph->getBigMDelay(), graph->getBigMJitter()), cshpGraph[k]);
    
    for(int i = 0; i < n; i++) 
        if (!graph->removed[i]) {
            add_vertex(i, graphEdmonds);
            add_vertex(i, heuristicGraph);
        }

    int countEdges = 0;
    for (int i = 0; i < n; i++)
        if(!graph->removed[i]) {
            for (auto arc : graph->arcs[i]) { 
                add_edge(i, arc->getD(), 0, graphEdmonds);
                add_edge(i, arc->getD(), arc->getDelay(), heuristicGraph);
                for (auto k : graph->DuS)
                    if (!graph->removedF[i][arc->getD()][k])
                        add_edge(i, arc->getD(), SPPRC_Graph_Arc(countEdges++, 0.0, arc->getDelay(), arc->getJitter()), cshpGraph[k]);
            }
        }
    cout << "Grafos criados" << endl;
    // getchar();
}

void Model::initialize() {
	int n = graph->getN();
	y = vector<vector<bool>>(n, vector<bool>(n));
    treeY = vector<vector<bool>>(n, vector<bool>(n));
    z = vector<bool>(n);
    f = vector<vector<vector<bool>>>(n, vector<vector<bool>>(n, vector<bool>(n)));
    used = vector<bool>(n);
}

bool Model::compareArcs(Arc_edmonds a, Arc_edmonds b) {
    return a.c < b.c;
}

void Model::edmonds() {
    Edge edge;
    int i, j, n = graph->getN(), root = graph->getRoot();
    double edgeWeight;
    y = vector<vector<bool >>(n, vector<bool>(n));
    vector<Edge> branching;
    arcsReturn = vector<Arc_edmonds>();

    Vertex roots[2];
    roots[0] = root;
	roots[1] = root;
    BOOST_FOREACH(Edge e, edges(graphEdmonds)) {
        arcsReturn.push_back(Arc_edmonds{int(e.m_source), int(e.m_target), get(weightMap, e)});
    }

    sort(arcsReturn.begin(), arcsReturn.end(), compareArcs);
    branchingEdges.erase(branchingEdges.begin(), branchingEdges.end());

    for (int i = 0; i < graph->getNAfterRemoved()-1; i++) {
        objectiveValue += arcsReturn[i].c;
        y[arcsReturn[i].o][arcsReturn[i].d] = true;   
    }

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

void Model::constrainedShortestpath(int k) {
    vector<vector<graph_traits<SPPRCGraph>::edge_descriptor>> opt_solutions;
    vector<spp_spp_res_cont> pareto_opt;
    int indexMin = 0, i, j, l;
    double minSP, minPath;
    vector<pair<int, int>> firstPath = vector<pair<int, int>>();

    SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, cshpGraph[k])[k];
    vert_prop.con_1 = graph->getParamDelay();
    vert_prop.con_2 = graph->getParamJitter();

    r_c_shortest_paths(cshpGraph[k],
                     get(&SPPRC_Graph_Vert::num, cshpGraph[k]),
                     get(&SPPRC_Graph_Arc::num, cshpGraph[k]),
                     graph->getRoot(),
                     k,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont(0, 0, 0),
                     ref_spprc(),
                     dominance_spptw(),
                     allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                     default_r_c_shortest_paths_visitor());

    if (!pareto_opt.empty()) {
        minSP = pareto_opt[0].cost;
        indexMin = 0;
        for (i = 1; i < opt_solutions.size(); i++) {
            if (pareto_opt[i].cost < minSP) {
                minSP = pareto_opt[i].cost;
                indexMin = i;
            }
        }

        minPath = minSP;
        for (j = static_cast<int>(opt_solutions[indexMin].size())-1; j >= 0; --j) {
            i = source(opt_solutions[indexMin][j], cshpGraph[k]);
            l = target(opt_solutions[indexMin][j], cshpGraph[k]);
            firstPath.push_back(make_pair(i, l));
        }
    } else z[k] = true;

    vert_prop.con_1 = graph->getBigMDelay();
    vert_prop.con_2 = graph->getBigMJitter();    

    r_c_shortest_paths(cshpGraph[k],
                     get(&SPPRC_Graph_Vert::num, cshpGraph[k]),
                     get(&SPPRC_Graph_Arc::num, cshpGraph[k]),
                     graph->getRoot(),
                     k,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont(0, 0, 0),
                     ref_spprc(),
                     dominance_spptw(),
                     allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                     default_r_c_shortest_paths_visitor()); 

    minSP = pareto_opt[0].cost;
    indexMin = 0;

    for (i = 1; i < opt_solutions.size(); i++) {
        if (pareto_opt[i].cost < minSP) {
            minSP = pareto_opt[i].cost;
            indexMin = i;
        }
    }

    if (z[k] || minSP < minPath) {
        for (j = static_cast<int>(opt_solutions[indexMin].size())-1; j >= 0; --j) {
            i = source(opt_solutions[indexMin][j], cshpGraph[k]);
            l = target(opt_solutions[indexMin][j], cshpGraph[k]);
            f[i][l][k] = true;
        }
        objectiveValue += minSP;
    } else {
        for (auto p : firstPath)
            f[p.first][p.second][k] = true;
        objectiveValue += minPath;
    }
    
}

void Model::updateEdgeBranching(int i, int j, double weight) {
    Edge e;
    bool found;
    tie(e, found) = edge(i, j, graphEdmonds);
    if (found) boost::put(edge_weight_t(), graphEdmonds, e, weight);
}

void Model::updateEdgePath(int i, int j, int k, double weight) {
    edge_descriptor ed;
    bool found;
    tie(ed, found) = edge(i, j, cshpGraph[k]);

    if (found) {
        SPPRC_Graph_Arc &arc_prop = get(edge_bundle, cshpGraph[k])[ed];
	    arc_prop.cost = weight;
	}
}

void Model::insertPenalties(double penalty) {
    objectiveValue += penalty;
}

bool Model::isAcyclic() {
    int n = graph->getN();

    BoostGraph G;

    for (int i = 0; i < n; i++)
        if (!graph->removed[i]) {
            // cout << i << ", ";
            add_vertex(i, G);
        }
    // cout << endl;

    for (int i = 0; i < graph->getNAfterRemoved()-1; i++) {
        // cout << arcsReturn[i].o << ", " << arcsReturn[i].d << endl;
        add_edge(arcsReturn[i].o, arcsReturn[i].d, 1, G);
    }

    // getchar();
    
    vector<int> component(n), discover_time(n);
    vector<default_color_type> color(n);
    vector<Vertex> rootG(n);

    int num = strong_components(G, make_iterator_property_map(component.begin(), get(vertex_index, G)),
                                root_map(make_iterator_property_map(rootG.begin(), get(vertex_index, G))).
                                        color_map(make_iterator_property_map(color.begin(), get(vertex_index, G))).
                                        discover_time_map(
                                        make_iterator_property_map(discover_time.begin(), get(vertex_index, G))));

    if (num != graph->getNAfterRemoved()) {
        cout << "Cycle: " << num << endl;
        return false;
    }
    return true;

}

int Model::heuristic() {
    int n = graph->getN();
    property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, heuristicGraph);
    vector<Vertex> predecessors = vector<Vertex>(n);
    vector<int> distance = vector<int>(n);
    vector<bool> notAttended = vector<bool>(n);
    vector<int> jitterDistance = vector<int>(n);

    dijkstra_shortest_paths(heuristicGraph, graph->getRoot(), predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, heuristicGraph))).distance_map(
            make_iterator_property_map(distance.begin(), get(vertex_index, heuristicGraph))));

    for (auto k : graph->terminals) 
        if (distance[k] > graph->getParamDelay()) 
            notAttended[k] = true;

    int actual, jitter, count = 0;
    for (auto t : graph->terminals) {
        actual = t, jitter = 0;
        while (actual != graph->getRoot()) {
            jitter += graph->getJitter(predecessors[actual], actual);
            actual = predecessors[actual];
            // cout << actual << " - " << graph->getRoot() << endl;
            // getchar();
        }
        if (jitter > graph->getParamJitter())
            notAttended[t] = true;
        jitterDistance[t] = jitter;
    }

    for (auto k : graph->terminals) {
        for (auto l : graph->terminals) {
            if (l != k && !notAttended[k] && !notAttended[l]) {
                if (distance[k] - distance[l] > graph->getParamVariation()) {
                    notAttended[distance[k] > distance[l] ? k : l] = true;
                    // notAttended[l] = true;
                }
            }
        }
    }

    for (auto t : graph->terminals) 
        if (notAttended[t])
            count++;  
    cout << "Heuristic: " << count << endl;  
    return count;
}