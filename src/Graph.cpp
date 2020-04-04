//
// Created by carlos on 27/05/19.
//

#include "Graph.h"
/*
// data structures for shortest path problem with resource constraint
// ResourceContainer model
struct spp_spp_res_cont_prep {
    spp_spp_res_cont(int c = 0, int r = 0) : cost(c), res(r) {}

    spp_spp_res_cont &operator=(const spp_spp_res_cont &other) {
        if (this == &other) return *this;
        this->~spp_spp_res_cont();
        new(this) spp_spp_res_cont(other);
        return *this;
    }

    int cost;
    int res;
};

bool operator==(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    return (res_cont_1.cost == res_cont_2.cost && res_cont_1.res == res_cont_2.res);
}

bool operator<(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    if (res_cont_1.cost > res_cont_2.cost) return false;
    if (res_cont_1.cost == res_cont_2.cost) return res_cont_1.res < res_cont_2.res;
    return true;
}

// ResourceExtensionFunction model
class ref_spprc_prep {
public:
    inline bool operator()(const SPPRCGraph &g, spp_spp_res_cont &new_cont, const spp_spp_res_cont &old_cont,
                           graph_traits<SPPRCGraph>::edge_descriptor ed) const {

        const SPPRC_Graph_Arc &arc_prop = get(edge_bundle, g)[ed];
        const SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.cost = old_cont.cost + arc_prop.cost;
        int &i_res = new_cont.res;
        i_res = old_cont.res + arc_prop.res;
        return i_res <= vert_prop.con;
    }
};

// DominanceFunction model
class dominance_spptw_prep {
public:
    inline bool operator()(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) const {
        return res_cont_1.cost <= res_cont_2.cost && res_cont_1.res <= res_cont_2.res;
    }
};
// end data structures for shortest path problem with time windows (spptw)
*/
Graph::Graph(string instance, string param, string outputName) {
    int u, v;
    double delay, jitter, bandwidth, ldp, paramDelayToken, paramJitterToken, paramVariationToken, paramBandwidthToken;
    int delayInt, jitterInt;
    string token;
    ifstream fileGraph, fileParam;
    ofstream output;

    output.open(outputName);

    fileParam.open(param, fstream::in);

    while (!fileParam.eof()) {
        fileParam >> token;
        if (token == "Delay") {
            fileParam >> token;
            if (token == "variation") {
                fileParam >> token >> paramVariationToken;
                Graph::paramVariation = int(1e5 * paramVariationToken);
            } else {
                fileParam >> paramDelayToken;
                Graph::paramDelay = int(1e5 * paramDelayToken);
            }
        }
        if (token == "Jitter") {
            fileParam >> token >> paramJitterToken;
            Graph::paramJitter = int(1e6 * paramJitterToken);
        }
        if (token == "Bandwidth") {
            fileParam >> token >> paramBandwidthToken;
            Graph::paramBandwidth = int(paramBandwidthToken);
        }
    }

    fileGraph.open(instance, fstream::in);

    while (!fileGraph.eof()) {
        fileGraph >> token;
        if (token == "Nodes") {
            fileGraph >> n;
            output << n << "\n";
            preProcessing = BoostGraph(n);
            arcs = vector<vector<Arc *>>(n, vector<Arc *>());
            removed = vector<bool>(n);
            removedF = vector<vector<vector<bool >>>(n, vector<vector<bool >>(n, vector<bool>(n)));

        }

        if (token == "Edges") {
            fileGraph >> m;
            output << m << "\n";
        }

        if (token == "E") {
            fileGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
            if (bandwidth >= paramBandwidth) {
                delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter), --u, --v;
                Arc *arc = new Arc(u, v, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                Arc *arcRev = new Arc(v, u, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                delayVector.push_back(delayInt), jitterVector.push_back(jitterInt);
                arcs[u].push_back(arc), arcs[v].push_back(arcRev);
                add_edge(u, v, delayInt, preProcessing), add_edge(v, u, delayInt, preProcessing);
            }
        }
        if (token == "Root") fileGraph >> root, root--;
        if (token == "T") fileGraph >> u, u--, terminals.push_back(u), DuS.push_back(u);
    }

    property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, preProcessing);
    predecessors = vector<VertexDescriptor>(n);
    distance = vector<int>(n);

    dijkstra_shortest_paths(preProcessing, root, predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, preProcessing))).distance_map(
            make_iterator_property_map(distance.begin(), get(vertex_index, preProcessing))));

    cntRemoved = n;
    for (int i = 0; i < n; i++) {
        // cout << i+1 << " - " << distance[i] << endl;
        removed[i] = distance[i] >= numeric_limits<int>::max();
        if (removed[i]) cntRemoved--;
    }
    output << cntRemoved << "\n";
    // getchar();

    bool isTerminal;
    for (int i = 0; i < n; ++i) {
        isTerminal = false;
        if (i != root) {
            for (auto t : terminals) {
                if (i == t) {
                    isTerminal = true;
                    break;
                }
            }
            if (!isTerminal) nonTerminals.push_back(i), DuS.push_back(i);
        }
    }

    sort(delayVector.begin(), delayVector.end(), greater<int>());
    sort(jitterVector.begin(), jitterVector.end(), greater<int>());

    for (int i = 0; i < cntRemoved - 1; i++)
        bigMDelay += delayVector[i], bigMJitter += jitterVector[i];

    output.close();
    cout << "Load graph successfully" << endl;
}
/*
void Graph::SAE() {
    int i, u, j, minSP, countEdges = 0;
    vector<int> jitterFromCShp = vector<int>(n), delayFromCShp = vector<int>(n), minSPVec = vector<int>(n);
    vector<int> distanceJitter = vector<int>(n);
    vector<vector<int>> CSHP = vector<vector<int>>(n, vector<int>(n));
    BoostGraph graphJitterSP = BoostGraph(n);
    SPPRCGraphPrep graphDelay, graphJitter;

    for (i = 0; i < n; i++)
        for (auto arc : arcs[i]) 
                add_edge(i, arc->getD(), arc->getJitter(), graphJitterSP);
                

    for (i = 0; i < n; i++) {
        distanceJitter = vector<int>(n);

        dijkstra_shortest_paths(graphJitterSP, i, predecessor_map(
                make_iterator_property_map(predecessors.begin(), get(vertex_index, graphJitterSP))).distance_map(
                make_iterator_property_map(distanceJitter.begin(), get(vertex_index, graphJitterSP))));

        minSP = paramJitter;
        for (auto t : terminals)
            if (t != i && distanceJitter[t] < minSP)
                minSP = distanceJitter[t];
        minSPVec[i] = minSP;

        add_vertex(SPPRC_Graph_Vert_Prep(i, paramJitter), graphDelay);
        add_vertex(SPPRC_Graph_Vert_Prep(i, paramDelay), graphJitter);
    }

    for (u = 0; u < n; ++u) {
        for (auto arc : arcs[u]) {
            j = arc->getD();
            add_edge(u, j, SPPRC_Graph_Arc_Prep(countEdges, arc->getDelay(), arc->getJitter()), graphDelay);
            add_edge(u, j, SPPRC_Graph_Arc_Prep(countEdges++, arc->getJitter(), arc->getDelay()), graphJitter);  
        }
    } 

    // Calculation of Constrained Shortests Paths
    vector<vector<graph_traits<SPPRCGraph>::edge_descriptor>> opt_solutions;
    vector<spp_spp_res_cont> pareto_opt;

    // CSP root -> S (NonTerminals)
    for (auto j : nonTerminals) {
        // CSHP delay
        if (!removed[j]) {
            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            vert_prop.con = paramJitter - minSPVec[j];

            r_c_shortest_paths(graphDelay,
                     get(&SPPRC_Graph_Vert_Prep::num, graphDelay),
                     get(&SPPRC_Graph_Arc_Prep::num, graphDelay),
                     root,
                     j,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont_prep(0, 0),
                     ref_spprc_prep(),
                     dominance_spptw_prep(),
                     allocator<r_c_shortest_paths_label<SPPRCGraphPrep, spp_spp_res_cont_prep >>(),
                     default_r_c_shortest_paths_visitor());
            if (pareto_opt.empty()) delayFromCShp[j] = paramDelay;
            else {
                minSP = pareto_opt[0].cost;
                for (auto p : pareto_opt) 
                    if (p.cost < minSP)
                        minSP = p.cost;
                delayFromCShp[j] = minSP;
            }
            vert_prop.con = paramJitter;

            SPPRC_Graph_Vert &vert_prop_d = get(vertex_bundle, graphJitter)[j];
            vert_prop_d.con = paramDelay;

            // CSHP jitter
            r_c_shortest_paths(graphJitter,
                     get(&SPPRC_Graph_Vert::num, graphJitter),
                     get(&SPPRC_Graph_Arc::num, graphJitter),
                     root,
                     j,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont(0, 0),
                     ref_spprc(),
                     dominance_spptw(),
                     allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                     default_r_c_shortest_paths_visitor());
            
            if (pareto_opt.empty()) jitterFromCShp[j] = paramJitter;
            else {
                minSP = pareto_opt[0].cost;
                for (auto p : pareto_opt) 
                    if (p.cost < minSP)
                        minSP = p.cost;
                jitterFromCShp[j] = minSP;
            }          
        }
    }

    // CSHP root -> j (terminals)
    for (auto j : terminals) {
        if (!removed[j]) {
            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            vert_prop.con = paramJitter;

            r_c_shortest_paths(graphDelay,
                     get(&SPPRC_Graph_Vert::num, graphDelay),
                     get(&SPPRC_Graph_Arc::num, graphDelay),
                     root,
                     j,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont(0, 0),
                     ref_spprc(),
                     dominance_spptw(),
                     allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                     default_r_c_shortest_paths_visitor());
            if (pareto_opt.empty()) delayFromCShp[j] = paramDelay;
            else {
                minSP = pareto_opt[0].cost;
                for (auto p : pareto_opt) 
                    if (p.cost < minSP)
                        minSP = p.cost;
                delayFromCShp[j] = minSP;
            }

            SPPRC_Graph_Vert &vert_prop_d = get(vertex_bundle, graphJitter)[j];
            vert_prop_d.con = paramDelay;
            // CSHP jitter
            r_c_shortest_paths(graphJitter,
                     get(&SPPRC_Graph_Vert::num, graphJitter),
                     get(&SPPRC_Graph_Arc::num, graphJitter),
                     root,
                     j,
                     opt_solutions,
                     pareto_opt,
                     spp_spp_res_cont(0, 0),
                     ref_spprc(),
                     dominance_spptw(),
                     allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                     default_r_c_shortest_paths_visitor());
            
            if (pareto_opt.empty()) jitterFromCShp[j] = paramJitter;
            else {
                minSP = pareto_opt[0].cost;
                for (auto p : pareto_opt) 
                    if (p.cost < minSP)
                        minSP = p.cost;
                jitterFromCShp[j] = minSP;
            }
        }
    }

    // CSHP: K -> J (NonTerminals) 
    for (auto j : DuS) {
        if (!removed[j]) {
            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            vert_prop.con = paramJitter - jitterFromCShp[j];

            for (auto k : terminals) {
                if (!removed[k] && k != j) {
                    r_c_shortest_paths(graphDelay,
                            get(&SPPRC_Graph_Vert::num, graphDelay),
                            get(&SPPRC_Graph_Arc::num, graphDelay),
                            k,
                            j,
                            opt_solutions,
                            pareto_opt,
                            spp_spp_res_cont(0, 0),
                            ref_spprc(),
                            dominance_spptw(),
                            allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                            default_r_c_shortest_paths_visitor());
                    if (pareto_opt.empty()) CSHP[k][j] = paramDelay;
                    else {
                        minSP = pareto_opt[0].cost;
                        for (auto p : pareto_opt) 
                            if (p.cost < minSP)
                                minSP = p.cost;
                        CSHP[k][j] = minSP;
                    }
                }
            }
            vert_prop.con = paramJitter;
        }
    }

    int cntRem = 0;
    for (auto i : DuS) {
        for (auto arc : arcs[i]) {
            j = arc->getD();
            for (auto k : terminals)
                if (k != j && k != i)
                    if (delayFromCShp[i] + arc->getDelay() + CSHP[k][j] > paramDelay) {
                        cntRem++;
                        removedF[i][j][k] = true;
                    }
        }
    }

    cout << "All paths are computed: " << cntRem << endl;
    // getchar();
}
*/
void Graph::showGraph() {
    cout << "Arcs" << endl;
    for (int o = 0; o < n; o++) {
        if (!removed[o]) {
            for (auto *arc : arcs[o])
                cout << arc->getO()+1 << " " << arc->getD()+1 << ": " << arc->getDelay()
                     << " " << arc->getJitter() << " " << arc->getBandwidth() << " " <<
                     arc->getEstimateLinkDuration() << endl;
        }
    }

    cout << "\n Param" << endl;
    cout << "Nodes: " << n << " Edges: " << m <<
         " CntTerminals: " << int(terminals.size()) << " Root: " << root << endl;

    cout << "Delay: " << paramDelay << " Jitter: " << paramJitter <<
         " DelayVari.: " << paramVariation << " Bandwidth: " << paramBandwidth << endl;

    cout << "Terminals" << endl;
    for (int i : terminals) {
        cout << "T: " << i+1 << " ";
    }
    cout << "\nNonTerminals" << endl;
    for (int i : nonTerminals) {
        cout << "NT: " << i+1 << " ";
    }
}

int Graph::getBigMDelay() {
    return bigMDelay;
}

int Graph::getBigMJitter() {
    return bigMJitter;
}

int Graph::getShpTerminal(int k) {
    return distance[k];
}

int Graph::getN() const {
    return n;
}

void Graph::setN(int n) {
    Graph::n = n;
}

int Graph::getM() const {
    return m;
}

void Graph::setM(int m) {
    Graph::m = m;
}

int Graph::getParamDelay() const {
    return paramDelay;
}

void Graph::setParamDelay(int paramDelay) {
    Graph::paramDelay = paramDelay;
}

int Graph::getParamJitter() const {
    return paramJitter;
}

void Graph::setParamJitter(int paramJitter) {
    Graph::paramJitter = paramJitter;
}

int Graph::getParamVariation() const {
    return paramVariation;
}

int Graph::getNAfterRemoved() {
    return cntRemoved;
}

void Graph::setParamVariation(int paramVariation) {
    Graph::paramVariation = paramVariation;
}

int Graph::getParamBandwidth() const {
    return paramBandwidth;
}

void Graph::setParamBandwidth(int paramBandwidth) {
    Graph::paramBandwidth = paramBandwidth;
}

int Graph::getRoot() const {
    return root;
}

void Graph::setRoot(int root) {
    Graph::root = root;
}

