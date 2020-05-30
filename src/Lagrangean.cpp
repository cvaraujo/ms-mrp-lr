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
    LB = 0, UB = model->heuristic(), iter = 0;
    // LB = 0, UB = int(graph->terminals.size()), iter = 0;
    cout << "Load Model" << endl;
}

void Lagrangean::getGradientTerminals(vector<vector<double>> &gradientVar) {
    int i, j, bigMK, bigML, n = graph->getN();
    for (int k: graph->terminals) {
        for (int l : graph->terminals) {
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                gradientVar[k][l] = -graph->getParamVariation() - model->z[k] * bigMK - model->z[l] * bigML;
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
    for (auto k : graph->terminals) {
        if (!graph->removed[k]) {
            for (i = 0; i < n; i++) {
                if (!graph->removed[i]) {
                    for (auto *arc : graph->arcs[i]) {
                        j = arc->getD();
                        gradientRel[i][j][k] = int(model->f[i][j][k]) - int(model->y[i][j]);
                        
                        if (gradientRel[i][j][k] > 0) {
                            feasible = false;
                        }
                    }
                }
            }
        }
    }
}

double Lagrangean::getNormRelation(vector<vector<vector<double>>> &gradient) {
    double sum = 0;
    int i, j, n = graph->getN();
    for (auto k : graph->terminals) {
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
        if (!graph->removed[i]) {
            for (auto *arc : graph->arcs[i]) {
                j = arc->getD();
                sumMultipliers = 0;
                for (auto k : graph->terminals) 
                    sumMultipliers -= multipliersRel[i][j][k];
                model->updateEdgeBranching(i, j, sumMultipliers);
            }
        }
    }
}

void Lagrangean::updatePathCosts(int k) {
    double edgeCost = 0, edgeCostAux;
    int i, j, n = graph->getN();

    for (i = 0; i < n; i++) {
        for (auto *arc : graph->arcs[i]) {
            j = arc->getD(), edgeCostAux = 0;
            if (!graph->removedF[i][j][k]) {
                for (auto l : graph->terminals)
                    if (!graph->removedF[i][j][l] && l != k) 
                        edgeCostAux += (multipliersVar[k][l] - multipliersVar[l][k]);
                
                edgeCost = multipliersRel[i][j][k] + (arc->getDelay() * edgeCostAux);
                model->updateEdgePath(i, j, k, edgeCost); 
            }               
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
                if (k != l) coefZ += multipliersVar[k][l] * bigMK + multipliersVar[l][k] * bigML;
            }
            
            if (1 - coefZ < 0 || model->noPath[k]) {
                model->z[k] = true;
            }
        }
    }

    for (int k: graph->terminals) {
        if (model->z[k]) {
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
    model->insertPenalties(penalty);
}

bool Lagrangean::solveModel() {
    model->objectiveValue = 0;
    model->initialize();
    
    for (auto k : graph->terminals) {
        if (!model->noPath[k]) {
            updatePathCosts(k);
            model->constrainedShortestpath(k);
        } else model->z[k] = true;
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
    if (feasible) return true;//model->isAcyclic();
    feasible = true;
    return false;
}

double Lagrangean::solve() {
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    int n = graph->getN(), originalObjectiveFunction, heuristicObj;
    iterBlb = 0, iterBub = 0, progress = 0;
    double thetaVar, normVar, thetaRel, normRel, objectiveFunctionPPL;

    vector<vector<double>> gradientVar;
    vector<vector<vector<double>>> gradientRel;

    gradientVar = vector<vector<double >>(n, vector<double>(n));
    gradientRel = vector<vector<vector<double>>>(graph->getN(), vector<vector<double>>(graph->getN(), vector<double>(graph->getN())));

    multipliersVar = vector<vector<double >>(n, vector<double>(n));
    multipliersRel = vector<vector<vector<double >>>(graph->getN(), vector<vector<double>>(graph->getN(), vector<double>(graph->getN())));

    // for (int i = 0; i < n; i++) {
        // cout << i << model->noPath[i] << endl;
    // }

    // getchar();

    while (iter < maxIter && endTime < time) {
        if (solveModel()) {
            getGradientTerminals(gradientVar);
            getGradientRelation(gradientRel);

            objectiveFunctionPPL = model->objectiveValue;
            // cout << "PPL: " << objectiveFunctionPPL << endl;

            if (objectiveFunctionPPL > LB) {
                LB = ceil(objectiveFunctionPPL), progress = 0;
                // LB = objectiveFunctionPPL, progress = 0;
                iterBlb = iter;
                if ((UB - LB) / UB <= 0.0001) return UB;
            } else { 
                progress++;
                if (progress >= B) {
                    lambda /= 2;
                    progress = 0;
                }
            }

            originalObjectiveFunction = originalObjectiveValue();
            heuristicObj = heuristic();
            // cout << "Original: " << originalObjectiveFunction << ", heuristic: " << heuristicObj << endl;

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
                    if (k != l) 
                        multipliersVar[k][l] = max(0.0, multipliersVar[k][l] + gradientVar[k][l] * thetaVar);

            for (auto k : graph->terminals)
                for (int i = 0; i < n; i++)
                    for (auto *arc : graph->arcs[i]) 
                        if (!graph->removedF[i][arc->getD()][k])
                            multipliersRel[i][arc->getD()][k] = max(0.0, multipliersRel[i][arc->getD()][k] + gradientRel[i][arc->getD()][k] * thetaRel);
            
            // cout << "(Feasible) Upper Bound = " << UB << ", (Relaxed) Lower Bound = " << LB << endl;

            iter++;
            end = chrono::steady_clock::now();
            endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
            // getchar();

        }
    }
    return 0;
}

int Lagrangean::heuristic() {
    // cout << "heuristic " << endl;
    int obj = 0, i, j, n = graph->getN(), bestCand;
    vector<int> delayPaths = vector<int>(n), jitterPaths = vector<int>(n),
                predecessors = vector<int>(n), delayPathAux = vector<int>(n),
                jitterPathAux = vector<int>(n), predecessorsAux = vector<int>(n);
    vector<bool> notAttended = vector<bool>(n), notAttendedAux = vector<bool>(n);
    vector<int> changed = vector<int>();

    predecessors[graph->getRoot()] = graph->getRoot();
    for (auto arc : model->branchingEdges) {
        predecessors[arc->getD()] = arc->getO();
    }

    random_shuffle(graph->terminals.begin(), graph->terminals.end());

    // cout << "Computing the paths" << endl;
    int actual, jitter, delay, count = 0;
    for (auto k : graph->DuS) {
        if (!graph->removed[k]) {
            actual = k, jitter = 0, delay = 0;
            while (actual != graph->getRoot()) {
                jitter += graph->getJitter(predecessors[actual], actual),
                delay += graph->getDelay(predecessors[actual], actual);
                actual = predecessors[actual];
            }
            delayPaths[k] = delay, jitterPaths[k] = jitter;
        }
    }

    obj = 0;
    for (auto k : graph->terminals) 
        if (delayPaths[k] > graph->getParamDelay() || jitterPaths[k] > graph->getParamJitter())
            notAttended[k] = true, obj++;

    // cout << "Before LS: " << obj << endl;

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
    if (selected == -1) return graph->terminals.size();

    // get the values of this path
    for (auto k : graph->terminals) {
        if (k != selected) {
            bestCand = -1, delay = graph->getParamDelay()+1, jitter = graph->getParamJitter()+1; 
            
            if (delayPaths[k] < (delayPaths[selected] - graph->getParamVariation()) || delayPaths[k] > (delayPaths[selected] + graph->getParamVariation())) {
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
                    notAttended[k] = false;
                    canMove = true;
                    losts = 0;
                    for (i = 0; i < n; i++) {
                        delayPathAux[i] = delayPaths[i], jitterPathAux[i] = jitterPaths[i], 
                        predecessorsAux[i] = predecessors[i], notAttendedAux[i] = notAttended[i];
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
                        for (int j : graph->DuS) {
                            if (!graph->removed[j] && j != actual && predecessorsAux[j] == actual) {
                                delayPathAux[j] += diffDelay, jitterPathAux[j] += diffJitter;
                                if (delayPathAux[j] > graph->getParamDelay() || jitterPathAux[j] > graph->getParamJitter()) {
                                    notAttendedAux[j] = true;
                                    losts++;
                                } else notAttendedAux[j] = false;

                                if (losts >= 2 || (j == selected && notAttendedAux[j])) {
                                    canMove = false;
                                    changed.erase(changed.begin(), changed.end());
                                    break;
                                } else changed.push_back(j);
                            }
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

    int leqSelec = 0, geqSelec = 0; 
    obj = 0;
    for (auto k : graph->terminals)
        if (k != selected) {
            if (notAttended[k] || delayPaths[k] < (delayPaths[selected] - graph->getParamVariation()) ||
                delayPaths[k] > (delayPaths[selected] + graph->getParamVariation()) ||
                delayPaths[k] > graph->getParamDelay() || jitterPaths[k] > graph->getParamJitter())
                obj++, notAttended[k] = true;
            
            if (!notAttended[k]) {
                if (delayPaths[k] <= delayPaths[selected]) leqSelec++;
                else geqSelec++;
            }
        } 

    // cout << leqSelec << " - " << geqSelec << endl;
    // getchar();
    for (auto k : graph->terminals) {
        if (k != selected && !notAttended[k]) {
            for (auto l : graph->terminals) {
                if (k != l && l != selected && !notAttended[l]) {
                    if (delayPaths[k] - delayPaths[l] > graph->getParamVariation()) {
                        if (leqSelec > geqSelec) {
                            if (delayPaths[k] > delayPaths[selected]) {
                                notAttended[k] = true;
                                obj++;
                                break;
                            } else notAttended[l] = true;
                        } else {
                            if (delayPaths[k] < delayPaths[selected]) {
                                notAttended[k] = true;
                                obj++;
                                break;
                            } else notAttended[l] = true;
                        }
                        obj++;
                    }
                }
            }
        }
    }

    // cout << "For Heuristic = " << obj << endl;
    // getchar();
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