#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <fstream>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>	Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Plane_3  Plane_3;

int main(int nargc, char *argv[]) {
	std::string s1filename, s2filename, soutfilename;
	Polyhedron P,Q,R;

	// Read in the two operand surfaces	
	if (nargc == 1) {
		std::cout << " Enter filename for E (as in I = E - F): ";
		std::cin >> s1filename;
		std::cout << " Enter filename for F (as in I = E - F): ";
		std::cin >> s2filename;
		std::cout << " Enter output file name I (as in I = E - F): ";
		std::cin >> soutfilename;
	}
	else if (nargc == 4) {
		s1filename = argv[1];
		s2filename = argv[2];
		soutfilename = argv[3];
	}
	else {
		std::cerr << " Wrong number of input arguments!" << std::endl;
	}

	std::ifstream infile1(s1filename.c_str());
	
	if (!infile1) {
		std::cerr << " Can not open the file: " << s1filename;
		exit(-1);
	}
		
	std::ifstream infile2(s2filename.c_str());
	
	if (!infile2) {
		std::cerr << " Can not open the file: " << s2filename;
		exit(-1);
	}
		
	infile1 >> P;
	infile2 >> Q;

	// Convert to Nef polyhedron
	Nef_polyhedron E(P);
	Nef_polyhedron F(Q);
	Nef_polyhedron G = Nef_polyhedron::COMPLETE;

	infile1.close();
	infile1.close();

	if (!P.is_closed()) {
		std::cout << "Stereo plane is not closed!" << std::endl;
	}

	// Perform subtraction boolean operation
	Nef_polyhedron I = E - F;

	if (I.is_simple()) {
		std::cout << " Diff is simple." << std::endl;
		I.convert_to_polyhedron(R);
		std::ofstream out2(soutfilename.c_str());
		out2 << R;
	}
	else {
		std::cout << " Diff is NOT simple." << std::endl;
		std::cout << " The resulting surface might not be sutiable for consequent meshing operation!" << std::endl << std::endl;
		std::cout << " Intersecting it with hyper-cube..." << std::endl;
		I = I * G;
		if (I.is_simple()) {
			std::cout << " I is simple." << std::endl;
		}
		I.convert_to_polyhedron(R);
		std::ofstream out2(soutfilename.c_str());
		out2 << R;
	}
	
	return 0;
}


