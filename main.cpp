#include <iostream>
#include <sstream>
#include <Graph.h>
#include <Lagrangean.h>

int main(int argc, const char *argv[]) {
	if (argc > 3) {
		// mkdir("results", 0777);
	    auto *graph = new Graph(argv[1], argv[2], argv[3]);
	    stringstream lambda(argv[4]);
	    stringstream maxIter(argv[5]);
	    stringstream B(argv[6]);
	    stringstream time(argv[7]);

	    double l;
	    int mIter, b, t;
	    lambda >> l;
	    maxIter >> mIter;
	    B >> b;
	    time >> t;

	    // cout << l << " - " << mIter << " - " << b << endl;
	    auto *lagrangean = new Lagrangean(graph, l, mIter, b, t);
	    lagrangean->solve();

	    lagrangean->showSolution(argv[3]);



	} else cout << "Wrong fields" << endl;
	
    return 0;
}