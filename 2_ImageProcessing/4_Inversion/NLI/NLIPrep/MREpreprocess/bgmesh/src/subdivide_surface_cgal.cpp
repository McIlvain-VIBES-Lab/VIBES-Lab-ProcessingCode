//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Subdivision_method_3.h>
#include <iostream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>

//typedef CGAL::Simple_cartesian<double>	   Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>		   Polyhedron;

using namespace std;
using namespace CGAL;

int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "Usage: refine_mesh d infilename outfilename" << endl;
		cout << "		d: the depth of the subdivision (0 < d < 10)" << endl;
		cout << "		filename: the input mesh (.off)" << endl;
		return 0;
	}
	
	
	int d = argv[1][0] - '0';
	
	std::ifstream ifs(argv[2]);
	Polyhedron P;
	ifs >> P; // read the .off

	Subdivision_method_3::Sqrt3_subdivision(P,d);

	std::ofstream ofs(argv[3]);
	ofs << P; // write the .off

	return 0;
}
