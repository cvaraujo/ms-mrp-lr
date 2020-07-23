//
// Created by carlos on 06/03/19.
//

#include <chrono>
#include "BarrierMethod.h"

BarrierMethod::BarrierMethod(Graph *graph) {
    if (graph != nullptr) {
        this->graph = graph;
    } else exit(EXIT_FAILURE);
}

void BarrierMethod::initialize() {
    int o, d, n = graph->getN(), m = graph->getM();
    try {
        env.set("LogFile", "MS_mip.log");
        env.start();

        f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
        y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
        z = vector<GRBVar>(n);

        char name[20];
        for (o = 0; o < n; o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                sprintf(name, "y_%d_%d", o, d);
                y[o][d] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
                for (int k: graph->DuS) {
                    sprintf(name, "f_%d_%d_%d", o, d, k);
                    if (graph->removedF[o][d][k]) this->f[o][d][k] = model.addVar(0.0, 0.0, 0, GRB_CONTINUOUS, name); 
                    else this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
                }
            }
        }

        for (auto i : graph->terminals) {
            sprintf(name, "z_%d", i);
            z[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
        }
        model.update();
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
        cout << ex.getErrorCode() << endl;
        exit(EXIT_FAILURE);
    }
}

void BarrierMethod::initModel(){
    cout << "BarrierMethod Creation" << endl;
    objectiveFunction();
    rootFlow(), flowConservation(), terminalsFlow();
    maxArcs(), relFandY();
    limDelayAndJitter(), limVariation();
    nonTerminalsLeafs(), primeToTerminals();
}

void BarrierMethod::objectiveFunction() {
    GRBLinExpr objective;
    for (auto k : graph->terminals)
        objective += z[k];

    model.setObjective(objective, GRB_MINIMIZE);
    model.update();
    cout << "Objective Function" << endl;
}

void BarrierMethod::rootFlow() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->terminals) {
        GRBLinExpr flowExpr, rootExpr;
        for (o = 0; o < graph->getN(); o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                if (o == root) flowExpr += f[root][d][k];
                else if (d == root) rootExpr += f[o][root][k];
            }
        }
        model.addConstr((flowExpr - rootExpr) == 1, "root_flow_all_" + to_string(k));
        beginConstrRel++;
    }
    model.update();
    cout << "Flow on root node" << endl;
}

void BarrierMethod::flowConservation() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->DuS) {
        for (int j = 0; j < graph->getN(); j++) {
            if (j != root && j != k) {
                GRBLinExpr flowIn, flowOut;
                for (o = 0; o < graph->getN(); o++) {
                    for (auto *arc : graph->arcs[o]) {
                        d = arc->getD();
                        if (o == j) flowOut += f[j][d][k];
                        if (d == j) flowIn += f[o][j][k];
                    }
                }
                model.addConstr((flowIn - flowOut) == 0, "flow_conservation_" + to_string(j) + "_" + to_string(k));
                beginConstrRel++;
            }
        }
    }
    model.update();
    cout << "Flow conservation" << endl;
}

void BarrierMethod::terminalsFlow() {
    int o, d;
    for (auto k : graph->DuS) {
        GRBLinExpr flowIn, flowOut;
        for (o = 0; o < graph->getN(); o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                if (o == k) flowOut += f[k][d][k];
                if (d == k) flowIn += f[o][k][k];
            }
        }
        model.addConstr((flowOut - flowIn) == -1, "flow_on_terminals_" + to_string(k));
        beginConstrRel++;
    }
    model.update(); 
    cout << "Flow on terminals" << endl;
}

void BarrierMethod::maxArcs() {
    GRBLinExpr totalArcs;
    for (int o = 0; o < graph->getN(); o++)
        for (auto *arc : graph->arcs[o]) {
            totalArcs += y[arc->getO()][arc->getD()];
        }

    model.addConstr(totalArcs == (graph->getN() - 1), "maximum_of_arcs");
    beginConstrRel++;
    model.update(); 
    cout << "maximum of arcs in the tree" << endl;
}

void BarrierMethod::relFandY() {
    endConstrRel = beginConstrRel;
    for (auto k : graph->DuS)
        for (int i = 0; i < graph->getN(); i++)
            for (auto *arc : graph->arcs[i]) {
                model.addConstr(f[i][arc->getD()][k] <= y[i][arc->getD()], "relation_" + to_string(i) + "_" + to_string(arc->getD()) + "_" + to_string(k));
                endConstrRel++;
            }
    beginConstrVar = endConstrRel;
    model.update();
}

void BarrierMethod::limDelayAndJitter() {
    int o, d, paramDelay, paramJitter;
    for (auto k : graph->terminals) {
        GRBLinExpr limDelay, limJitter;
        for (o = 0; o < graph->getN(); o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                limDelay += arc->getDelay() * f[o][d][k];
                limJitter += arc->getJitter() * f[o][d][k];
            }
        }

        paramDelay = graph->getParamDelay(), paramJitter = graph->getParamJitter();
        model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]), "delay_limit_" + to_string(k));
        model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]), "jitter_limit_" + to_string(k));
        beginConstrVar+=2;
    }
    model.update();
    cout << "Delay and Jitter limits" << endl;
}

void BarrierMethod::limVariation() {
    endConstrVar = beginConstrVar;
    int o, d, bigMK, bigML;
    for (auto k : graph->terminals) {
        for (auto l : graph->terminals) {
            if (k != l) {
                GRBLinExpr delayVariation;
                for (o = 0; o < graph->getN(); o++) {
                    for (auto *arc : graph->arcs[o]) {
                        d = arc->getD();
                        delayVariation += arc->getDelay() * (f[o][d][k] - f[o][d][l]);
                    }
                }
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                
                model.addConstr(delayVariation <= graph->getParamVariation() + bigMK * z[k] + bigML * z[l],
                                "limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
                endConstrVar++;
            }
        }
    }
    model.update();
    cout << "Delay variation limits" << endl;
}

void BarrierMethod::nonTerminalsLeafs() {
    beginConstrLeaf = endConstrVar;
    endConstrLeaf = beginConstrLeaf;
    for (auto q : graph->DuS) {
        for (auto e : graph->DuS) {
            if (e != q) {
                model.addConstr(f[0][q][e] <= 0, "non_terminals_leafs_" + to_string(q) + "_" + to_string(e));
                endConstrLeaf++;
            }
        }
    }
    model.update();
    cout << "Non terminals leafs" << endl;
}

void BarrierMethod::primeToTerminals() {
    model.addConstr(y[graph->getRoot()][0] == 1);
    for (auto k : graph->terminals) 
        model.addConstr(z[k] >= f[0][k][k], "prime_to_terminals_" + to_string(k));

    model.update();
    cout << "S' to terminals" << endl;
}

void BarrierMethod::solve() {
    try {
        model.set("TimeLimit", "3600.0");
        model.set("Method", "2");
        model.set("Crossover", "0");
        model.update();
        model.write("linear_model.lp");
        model.optimize();
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
        exit(0);
    }
}

void BarrierMethod::barrierMethod() {
    int i, j, k, bigMK, bigML, n = graph->getN();
    initialize();
    initModel();
    solve();

    cout << "Obj LR: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    multipliersVar = vector<vector<double >>(n, vector<double>(n));
    multipliersLeaf = vector<vector<double >>(n, vector<double>(n));
    multipliersRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));

    for (int c = beginConstrRel; c < endConstrRel; c++) {
        auto constr = model.getConstr(c);
        string text = constr.get(GRB_StringAttr_ConstrName);
        vector<string> result;
        boost::split(result, text, boost::is_any_of("_"));
        stringstream ic(result[1]), jc(result[2]), kc(result[3]);
        ic >> i; jc >> j ; kc >> k;
        multipliersRel[i][j][k] = abs(constr.get(GRB_DoubleAttr_Pi));
    }

    for (int c = beginConstrVar; c < endConstrVar; c++) {
        auto constr = model.getConstr(c);
        string text = constr.get(GRB_StringAttr_ConstrName);
        vector<string> result;
        boost::split(result, text, boost::is_any_of("_"));
        stringstream ic(result[5]), kc(result[6]);
        ic >> i; kc >> k;
        multipliersVar[k][i] = abs(constr.get(GRB_DoubleAttr_Pi));
    }

    for (int c = beginConstrLeaf; c < endConstrLeaf; c++) {
        auto constr = model.getConstr(c);
        string text = constr.get(GRB_StringAttr_ConstrName);
        vector<string> result;
        boost::split(result, text, boost::is_any_of("_"));
        stringstream ic(result[3]), kc(result[4]);
        ic >> i; kc >> k;
        multipliersLeaf[i][k] = abs(constr.get(GRB_DoubleAttr_Pi));
    }

    cout << "Finished the update" << endl;
}
