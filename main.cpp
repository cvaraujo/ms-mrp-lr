#include <iostream>
#include <sstream>
#include <Graph.h>
#include <Lagrangean.h>

int main(int argc, const char *argv[]) {
  if (argc >= 12) {
    // to random decisions
    srand(time(NULL));

    // Create the graph and fill the parameters 
    auto *graph = new Graph(argv[2], argv[3], argv[4]);
	    
    // Preprocessing
    stringstream prepMve(argv[5]);
    stringstream prepSae(argv[6]);
    bool mve, sae;
    prepMve >> mve; prepSae >> sae;

    if (mve) graph->MVE(argv[4]);
    if (sae) graph->SAE(argv[4]);
    graph->finishPreprocessing(argv[4], mve, sae);
    // ./MSLagrangean relaxation graph.txt param.txt result.txt
    // MVE SAE barrierMethod heuristics 
    // stepSize maxIter B timeLimit 
	    
    stringstream rlSelected(argv[1]);
    stringstream barrier(argv[7]);
    stringstream heuristics(argv[8]);
    int rl;
    bool bm, useHeuristics;

    rlSelected >> rl;
    barrier >> bm;
    heuristics >> useHeuristics;
	    
    // Parameters of lagrangean relaxation
    stringstream lambda(argv[9]);
    stringstream maxIter(argv[10]);
    stringstream B(argv[11]);
    stringstream time(argv[12]);

    double l;
    int mIter, b, t;
    lambda >> l; maxIter >> mIter; B >> b; time >> t;

    auto *lagrangean = new Lagrangean(graph, rl, useHeuristics, bm, l, mIter, b, t);
    lagrangean->solve();
    lagrangean->showSolution(argv[4]);
  } else {
    cout << "Wrong number of parameters, try: \n ./MSLagrangean graph-file.txt param-file.txt result-file.txt " << endl;
	
  }
  return 0;
}
