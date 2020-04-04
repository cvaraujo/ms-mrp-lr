//
// Created by carlos on 31/05/19.
//

#include "Lagrangean.h"

Lagrangean::Lagrangean(Graph *graph, double lambda, int maxIter, int B, int time) {
    Lagrangean::graph = graph;
    Lagrangean::lambda = lambda;
    Lagrangean::maxIter = maxIter;
    Lagrangean::B = B;
    Lagrangean::time = time;
    model = new Model(graph);
    LB = 0, UB = int(graph->terminals.size()), iter = 0;
    cout << "Load Model" << endl;
}

void Lagrangean::getGradientTerminals(vector<vector<double>> &gradientVar) {
    int i, j, bigMK, bigML, n = graph->getN();
    for (int k: graph->terminals) {
        for (int l : graph->terminals) {
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                gradientVar[k][l] = -graph->getParamVariation() - (model->z[k] * bigMK + model->z[l] * bigML);
                for (i = 0; i < n; i++) {
                    for (auto *arc : graph->arcs[i]) {
                        j = arc->getD();
                        gradientVar[k][l] += arc->getDelay() * (int(model->f[i][j][k]) - int(model->f[i][j][l]));
                    }
                }
                if (gradientVar[k][l] > 0) feasible = false;
            }
        }
    }
}

void Lagrangean::getGradientRelation(vector<vector<vector<double>>> &gradientRel) {
    int i, j, n = graph->getN();
    for (auto k : graph->DuS) {
        for (i = 0; i < n; i++) {
            for (auto *arc : graph->arcs[i]) {
                j = arc->getD();
                gradientRel[i][j][k] = int(model->f[i][j][k]) - int(model->y[i][j]);
                if (gradientRel[i][j][k] > 0) feasible = false;
            }
        }
    }
}

double Lagrangean::getNormRelation(vector<vector<vector<double>>> &gradient) {
    double sum = 0;
    int i, j, n = graph->getN();
    for (auto k : graph->DuS) {
        for (i = 0; i < n; i++) {
            for (auto *arc : graph->arcs[i]) {
                j = arc->getD();
                sum += pow(gradient[i][j][k], 2);
            }
        }
    }
    return sqrt(sum);
}

double Lagrangean::getNormTerminals(vector<vector<double>> &gradient) {
    double sum = 0;
    for (int k : graph->terminals)
        for (int l : graph->terminals)
            if (k != l) sum += pow(gradient[k][l], 2);
    return sqrt(sum);
}

void Lagrangean::updateTreeCosts() {
    int i, j, n = graph->getN();
    double sumMultipliers;
    for (i = 0; i < n; i++) {
        for (auto *arc : graph->arcs[i]) {
            j = arc->getD();
            sumMultipliers = 0;
            for (auto k : graph->DuS) sumMultipliers -= multipliersRel[i][j][k];
            model->updateEdgeBranching(i, j, sumMultipliers);
        }
    }
}

void Lagrangean::updatePathCostsNT() {
    int i, j;
    double sumMultipliers;

    for(i = 0; i < graph->getN(); i++) {
        for (auto arc : graph->arcs[i]) {
            j = arc->getD();
            sumMultipliers = 0;
            for (auto q : graph->nonTerminals)
                sumMultipliers += multipliersRel[i][j][q];
            
            model->updateEdgePath(i, j, sumMultipliers, false);
        }               
    }
}

void Lagrangean::updatePathCosts(int k) {
    double edgeCost = 0, edgeCostAux;
    int i, j, n = graph->getN();

    for (i = 0; i < n; i++) {
        for (auto *arc : graph->arcs[i]) {
            j = arc->getD(), edgeCostAux = 0;
            for (auto l : graph->terminals) {
                if (l != k) edgeCostAux += (multipliersVar[k][l] - multipliersVar[l][k]);
            }
            edgeCost = multipliersRel[i][j][k] + (arc->getDelay() * edgeCostAux);
            model->updateEdgePath(i, j, edgeCost, true);                
        }
    }

}

void Lagrangean::updatePPL() {
    double penalty = 0, coefZ;
    int bigML, bigMK;

    for (int k : graph->terminals) {
        if (!model->z[k]) {
            coefZ = 0;
            for (int l : graph->terminals) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                if (k != l) {
                    coefZ += multipliersVar[k][l] * bigMK;
                    coefZ += multipliersVar[l][k] * bigML;
                }
            }
            if (1 - coefZ < 0) {
                model->z[k] = true;
            }
        }
    }

    for (int k: graph->terminals) {
        if (model->z[k]) {
            // cout << k << ", ";
            penalty += 1;
        }
        for (int l : graph->terminals) 
            if (k != l) {
                penalty -= (multipliersVar[k][l] * graph->getParamVariation());
                if (model->z[k]) {
                    bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                    bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                    penalty -= (multipliersVar[k][l] * bigMK) + (multipliersVar[l][k] * bigML);
                }
            } 
    }
    // cout << endl;
    // getchar();
    // cout << "Penalty: " << penalty << endl;
    model->insertPenalties(penalty);
}

bool Lagrangean::solveModel() {
    model->objectiveValue = 0;
    model->initialize();
    for (auto k : graph->terminals) {
        updatePathCostsNT();
        // cout << "NT" << endl;
        updatePathCosts(k);
        // cout << "T" << endl;
        // getchar();
        model->constrainedShortestpath(k);
        // cout << "PATH" << endl;

    }

    // cout << "CSHP: " << model->objectiveValue << endl;

    updateTreeCosts();
    model->edmonds();

    // cout << "tree: " << model->objectiveValue << endl;

    updatePPL();

    // cout << "Penalty: " << model->objectiveValue << endl;

    return true;
}

double Lagrangean::originalObjectiveValue() {
    int n = graph->getN(), j;
/*
    for (int i = 0; i < n; i++) {
        for (auto arc : graph->arcs[i]) {
            j = arc->getD();
            if(model->y[i][j]) {
                cout << i << " - " << j << endl;
            }
        }
    }
    cout << "---------------------" << endl;
    for (auto k : graph->terminals) {
        for (int i = 0; i < n; i++) {
            for (auto arc : graph->arcs[i]) {
                j = arc->getD();
                if(model->f[i][j][k]) {
                    cout << i << " - " << j << " - " << k << endl;
                }
            }
        }
    }
    cout << "---------------------" << endl;
*/
    int obj = 0;
    for (auto k : graph->terminals)
        if (model->z[k]) {
            // cout << k << endl;
            obj++;
        }
    return obj;
}

bool Lagrangean::isFeasible() {
    if (feasible) return true;
    feasible = true;
    return false;
}

double Lagrangean::solve() {
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    int n = graph->getN(), originalObjectiveFunction, heuristicObj;
    iterBlb = 0, iterBub = 0;
    double thetaVar, normVar, thetaRel, normRel, objectiveFunctionPPL;

    vector<vector<double>> gradientVar;
    vector<vector<vector<double>>> gradientRel;

    gradientVar = vector<vector<double >>(n, vector<double>(n));
    gradientRel = vector<vector<vector<double>>>(graph->getN(), vector<vector<double>>(graph->getN(), vector<double>(graph->getN())));

    multipliersVar = vector<vector<double >>(n, vector<double>(n));
    multipliersRel = vector<vector<vector<double >>>(graph->getN(), vector<vector<double>>(graph->getN(), vector<double>(graph->getN())));

    while (iter < maxIter && endTime < time) {
        if (solveModel()) {
            getGradientTerminals(gradientVar);
            getGradientRelation(gradientRel);

            objectiveFunctionPPL = model->objectiveValue;
            cout << "PPL: " << objectiveFunctionPPL << endl;

            if (objectiveFunctionPPL > LB) {
                LB = ceil(objectiveFunctionPPL), progress = 0;
                iterBlb = iter;
                if ((UB - LB) / UB <= 0.0001) return UB;
            } else { 
                progress++;
                if (progress == B) {
                    lambda /= 2;
                    progress == 0;
                }
            }

            originalObjectiveFunction = originalObjectiveValue();
            heuristicObj = heuristic();
            cout << "Original: " << originalObjectiveFunction << ", heuristic: " << heuristicObj << endl;

            if (isFeasible() && originalObjectiveFunction < UB) {
                UB = originalObjectiveFunction;
                iterBub = iter;
                if ((UB - LB) / UB <= 0.0001) return UB;
            } else if(heuristicObj < UB) {
                UB = heuristicObj;
                iterBub = iter;
                if ((UB - LB) / UB <= 0.0001) return UB;
            }

            normVar = getNormTerminals(gradientVar);
            normRel = getNormRelation(gradientRel);

            if (normVar == 0) thetaVar = 0;
            else thetaVar = lambda * ((UB - objectiveFunctionPPL) / pow(normVar, 2));

            if (normRel == 0) thetaRel = 0;
            else thetaRel = lambda * ((UB - objectiveFunctionPPL) / pow(normRel, 2));

            for (int k : graph->terminals) 
                for (int l : graph->terminals)
                    if (k != l) {
                        multipliersVar[k][l] = max(0.0, multipliersVar[k][l] + gradientVar[k][l] * thetaVar);
                    }
            
            for (auto k : graph->DuS)
                for (int i = 0; i < n; i++)
                    for (auto *arc : graph->arcs[i]) 
                        multipliersRel[i][arc->getD()][k] = max(0.0, multipliersRel[i][arc->getD()][k] + gradientRel[i][arc->getD()][k] * thetaRel);
            
            cout << "(Feasible) Upper Bound = " << UB << ", (Relaxed) Lower Bound = " << LB << endl;

            iter++;
            end = chrono::steady_clock::now();
            endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
            // getchar();

        }
    }
    return 0;
}

int Lagrangean::heuristic() {
    // cout << "heuristic ";
    int obj = 0, i, j, n = graph->getN();
    vector<vector<Arc*>> tree = vector<vector<Arc *>>(n, vector<Arc *>());
    vector<int> delayPaths = vector<int>(n);
    vector<int> jitterPaths = vector<int>(n);
    vector<bool> auxZ = vector<bool>(n);
    
    for (i = 0; i < n; i++) {
        for (auto arc : graph->arcs[i]) {
            if (model->treeY[i][arc->getD()]) {
                j = arc->getD();
                Arc *aux = new Arc(j, i, arc->getDelay(), arc->getJitter(), 0, 0);
                tree[j].push_back(aux);
            }
        }
    }
    
   for (auto k : graph->terminals) {
        i = tree[k][0]->getD();
        delayPaths[k] += tree[k][0]->getDelay();
        jitterPaths[k] += tree[k][0]->getJitter();
        while (i != graph->getRoot()) {
            delayPaths[k] += tree[i][0]->getDelay();
            jitterPaths[k] += tree[i][0]->getJitter();
            i = tree[i][0]->getD();
        }
    }

    for (auto k : graph->terminals) {
        if (delayPaths[k] > graph->getParamDelay() || jitterPaths[k] > graph->getParamJitter()) {
            auxZ[k] = true;
        }
    }

    int aux = 0;

    for (auto k : graph->terminals)
        if (auxZ[k]) 
            aux;

    // cout << "Only path = " << aux << endl;

    for (auto k : graph->terminals) 
        if (!auxZ[k]) 
            for (auto l : graph->terminals) 
                if (!auxZ[l]) 
                    if (k != l) 
                        if (delayPaths[k] - delayPaths[l] > graph->getParamVariation()) 
                            auxZ[k] = true;

    for (auto k : graph->terminals)
        if (auxZ[k]) 
            obj++;

    cout << "For Variation = " << obj << endl;
    // getchar();
    // cout << "finished" << endl;
    return obj;
}

void Lagrangean::showSolution(string outputName) {
    ofstream output;
    output.open(outputName);

    output << lambda << " " << maxIter << " " << B << " " << time << endl;

    output << UB << " " << LB << " " << 100 * ((UB - LB) / UB) << endl;

    output << iterBlb << " " << iterBub << " " << endTime << endl;

    output.close();
}