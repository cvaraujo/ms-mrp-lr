//
// Created by carlos on 31/05/19.
//

#include "../headers/Lagrangean.h"

Lagrangean::Lagrangean(Graph *graph, int relaxNum, bool heuristics, bool barrierMethod, double lambda, int maxIter, int B, int time) {
  Lagrangean::graph = graph;
  Lagrangean::lambda = lambda;
  Lagrangean::maxIter = maxIter;
  Lagrangean::B = B;
  Lagrangean::time = time;
  Lagrangean::relaxNum = relaxNum;
  Lagrangean::heuristics = heuristics;

  // Load the model: Graphs for paths and heuristics
  model = new Model(graph, relaxNum, heuristics);

  LB = 0, UB = graph->terminals.size(), iter = 0;
  if (heuristics) UB = model->initialHeuristic();
  firstUB = UB;
  
  int n = graph->getN();
  multipliersDelay = vector<double>(n);
  multipliersJitter = vector<double>(n);
  multipliersVar = vector<vector<double >>(n, vector<double>(n));
  multipliersLeaf = vector<vector<double >>(n, vector<double>(n));
  multipliersRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));
  freqLB = vector<vector<int >>(n, vector<int>(n));
  freqUB = vector<vector<int >>(n, vector<int>(n));

  if (barrierMethod) {
    auto start = chrono::steady_clock::now();
    BarrierMethod *bm  = new BarrierMethod(graph);
    bm->initModel();
    bm->solve();
    if (relaxNum == 1 || relaxNum == 3) bm->getMultipliersDelay(multipliersDelay);
    if (relaxNum <= 2) bm->getMultipliersJitter(multipliersJitter);

    bm->getMultipliersRelation(multipliersRel);
    bm->getMultipliersVariation(multipliersVar);
    bm->getMultipliersLeaf(multipliersLeaf);

    auto end  = chrono::steady_clock::now();
    bmTime = chrono::duration_cast<chrono::seconds>(end - start).count();
  }
  cout << "Lagrangean was initialized!" << endl;
}

void Lagrangean::getGradientDelay(vector<double> &gradientDelay){
  int i, j;
  for (auto k : graph->terminals) {
    gradientDelay[k] = -(graph->getParamDelay() + model->z[k] * graph->getBigMDelay()); 
    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	gradientDelay[k] += arc->getDelay() * model->f[i][j][k];
      }
    }
    if (gradientDelay[k] > 0) feasible = false;
  }
}

void Lagrangean::getGradientJitter(vector<double> &gradientJitter){
  int i, j;
  for (auto k : graph->terminals) {
    gradientJitter[k] = -(graph->getParamJitter() + model->z[k] * graph->getBigMJitter()); 
    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	gradientJitter[k] += arc->getJitter() * model->f[i][j][k];
      }
    }
    if (gradientJitter[k] > 0) feasible = false;
  }
}

void Lagrangean::getGradientVariation(vector<vector<double>> &gradientVar) {
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
	    gradientVar[k][l] += arc->getDelay() * (model->f[i][j][k] - model->f[i][j][l]);
	  }
	}
	if (gradientVar[k][l] > 0) feasible = false;
      }
    }
  }
}

void Lagrangean::getGradientLeaf(vector<vector<double>> &gradientLeaf) {
  for (auto q : graph->DuS) {
    for (auto e : graph->DuS) {
      if (q != e) {
	gradientLeaf[q][e] = int(model->f[0][q][e]);
	if (gradientLeaf[q][e] > 0) feasible = false;
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

double Lagrangean::getNormDelay(vector<double> &gradientDelay) {
  double sum = 0;
  for (auto k : graph->terminals)
    sum += pow(gradientDelay[k], 2);
  return sqrt(sum);
}

double Lagrangean::getNormJitter(vector<double> &gradientJitter) {
  double sum = 0;
  for (auto k : graph->terminals)
    sum += pow(gradientJitter[k], 2);
  return sqrt(sum);
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

double Lagrangean::getNormLeaf(vector<vector<double>> &gradient) {
  double sum = 0;
  for (int q : graph->DuS)
    for (int e : graph->DuS)
      if (q != e) sum += pow(gradient[q][e], 2);
  return sqrt(sum);
}

double Lagrangean::getNormVariation(vector<vector<double>> &gradient) {
  double sum = 0;
  for (int k : graph->terminals)
    for (int l : graph->terminals)
      if (k != l) sum += pow(gradient[k][l], 2);
  return sqrt(sum);
}

bool Lagrangean::isFeasible() {
  if (feasible) return model->isAcyclic();
  else feasible = true;
  return false;
}

double Lagrangean::solve() {
  // Time limit
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    
  int n = graph->getN(), originalObj, heuristicObj = int(graph->terminals.size());
  iterBlb = 0, iterBub = 0, progress = 0;
    
  double thetaVar, normVar, thetaRel, normRel,  thetaDelay, normDelay, thetaJitter, normJitter, lpObj;

  // Auxiliar vectors
  vector<double> gradientDelay, gradientJitter;
  vector<vector<double>> gradientVar;
  vector<vector<vector<double>>> gradientRel;

  gradientDelay = vector<double>(n);
  gradientJitter = vector<double>(n);
  gradientVar = vector<vector<double >>(n, vector<double>(n));
  gradientRel = vector<vector<vector<double>>>(n, vector<vector<double>>(n, vector<double>(n)));

  while (iter < maxIter && endTime < time) {
    if (model->solve(multipliersDelay, multipliersJitter, multipliersVar, multipliersRel)) {
      if (iter ==  0) firstLB = model->getObj();
      if (relaxNum == 1 || relaxNum == 3) getGradientDelay(gradientDelay);
      if (relaxNum <= 2) getGradientJitter(gradientJitter);

      getGradientVariation(gradientVar);
      getGradientRelation(gradientRel);

      lpObj = model->getObj();
      // Improvement of the lower bound?
      if (lpObj > LB) {
	LB = lpObj, progress = 0;
	iterBlb = iter;
      } else { 
	progress++;
	if (progress >= B) {
	  lambda /= 2;
	  progress = 0;
	}
      }
      originalObj = model->getOriginalObj();

      if (isFeasible() && originalObj < UB) {
	UB = originalObj, iterBub = iter;
	if ((UB - LB) / UB <= 0.0001) return UB;
      } 

      if (heuristics) {
	heuristicObj = model->getHeuristicObj();
	if(heuristicObj < UB) {
	  UB = heuristicObj, iterBub = iter;
	  if ((UB - LB) / UB <= 0.0001) return UB;
	}
      }

      for (int i = 0; i < graph->getN(); i++) {
	for (auto *arc : graph->arcs[i]) {
	  for (auto k : graph->DuS)
	    if (model->f[i][arc->getD()][k])
	      freqLB[i][arc->getD()] += 1;
	  if (model->treeY[i][arc->getD()])
	    freqUB[i][arc->getD()] += 1;
	}
      }
      //      cout << "PPL: " << lpObj << ", Original: " << originalObj << ", heuristic: " << heuristicObj << endl;
      // Step size 
      if (relaxNum == 1 || relaxNum == 3){
	normDelay = getNormDelay(gradientDelay);
	if (normDelay == 0) thetaDelay = 0;
	else thetaDelay = lambda * ((UB - lpObj) / pow(normDelay, 2));
      }

      if (relaxNum <= 2) {
	normJitter = getNormJitter(gradientJitter);
	if (normJitter == 0) thetaJitter = 0;
	else thetaJitter = lambda * ((UB - lpObj) / pow(normJitter, 2));
      }
      
      normVar = getNormVariation(gradientVar);
      normRel = getNormRelation(gradientRel);

      if (normVar == 0) thetaVar = 0;
      else thetaVar = lambda * ((UB - lpObj) / pow(normVar, 2));

      if (normRel == 0) thetaRel = 0;
      else thetaRel = lambda * ((UB - lpObj) / pow(normRel, 2));

      // Update the multipliers
      for (int k : graph->terminals) { 
	if (relaxNum == 1 || relaxNum == 3) multipliersDelay[k] = max(0.0, multipliersDelay[k] + gradientDelay[k] * thetaDelay);
	if (relaxNum <= 2) multipliersJitter[k] = max(0.0, multipliersJitter[k] + gradientJitter[k] * thetaJitter);
	for (int l : graph->terminals)
	  if (k != l) multipliersVar[k][l] = max(0.0, multipliersVar[k][l] + gradientVar[k][l] * thetaVar);
      }
      
      for (auto k : graph->DuS)
	for (int i = 0; i < n; i++)
	  for (auto *arc : graph->arcs[i]) 
	    multipliersRel[i][arc->getD()][k] = max(0.0, multipliersRel[i][arc->getD()][k] + (gradientRel[i][arc->getD()][k] * thetaRel));
    
      cout << "(Feasible) Upper Bound = " << UB << ", (Relaxed) Lower Bound = " << LB << endl;
      //getchar();
      iter++;
      end = chrono::steady_clock::now();
      endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    }
  }
  return 0;
}

void Lagrangean::showSolution(string outputName) {
  ofstream output;
  output.open(outputName, ofstream::app);

  output << "First LB: " << firstLB << "\nLB: " << LB << "\nIter. LB: " << iterBlb << endl;
  output << "First UB: " << firstUB << "\nUB: " << UB << "\nIter. UB: " << iterBub << endl;
    
  if (LB < 0) LB = 0;
  
  output << "gap: " << 100 * (double(UB - ceil(LB)) / double(UB)) << endl;
  output << "BM. Time: " << bmTime << "\nRuntime: " << endTime << endl;

  for (int i = 0; i < graph->getN(); i++)
    for (auto *arc : graph->arcs[i])
      output << "FL " << i << " " << arc->getD() << " " << freqLB[i][arc->getD()] << endl;

  for (int i = 0; i < graph->getN(); i++)
    for (auto *arc : graph->arcs[i])
      output << "FU " << i << " " << arc->getD() << " " << freqUB[i][arc->getD()] << endl;

  output.close();
}
