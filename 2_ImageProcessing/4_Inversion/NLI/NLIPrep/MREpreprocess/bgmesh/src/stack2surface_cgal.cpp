#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

// Usage:
// stack2surface_cgal infn.inr outfn.off [parameter.txt]
int main(int argc, char *argv[]) {
	Tr tr;			// 3D-Delaunay triangulation
	C2t3 c2t3 (tr);	// 2D-complex in 3D-Delaunay triangulation
	
	double isovalue;
	double scx;
	double scy;
	double scz;
	double bssr;
	double prec;
	double facet_angle;
	double facet_radius;
	double facet_distance;
	
	if (argc==4) {
		std::ifstream ifs(argv[3]);
		ifs >> isovalue;
		ifs >> scx;
		ifs >> scy;
		ifs >> scz;
		ifs >> bssr;
		ifs >> prec;
		ifs >> facet_angle;
		ifs >> facet_radius;
		ifs >> facet_distance;
	}
	else if (argc == 3){
		isovalue=0.9;
		scx=125; scy=125; scz=62;
		bssr=254*254*2.;
		prec=1e-8;
		facet_angle=30.;
		facet_radius=3.;
		facet_distance=3.;
	}
	else {
		std::cerr << " Need at least 3 input arguments\n";
		std::cerr << " Usage: stack2surface_cgal infn.inr outfn.off [parameter.txt]\n";
		exit(1);
	}
	  // the 'function' is a 3D gray level image
	std::string inrfn(argv[1]);
	Gray_level_image image(inrfn.c_str(), isovalue);
	

	GT::Point_3 bounding_sphere_center(scx, scy, scz);
	GT::FT bounding_sphere_squared_radius = bssr;
	GT::Sphere_3 bounding_sphere(bounding_sphere_center,
								   bounding_sphere_squared_radius);

	Surface_3 surface(image, bounding_sphere, prec);

	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(facet_angle,
													 facet_radius,
													 facet_distance);

	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
	std::string outfn(argv[2]);
	std::ofstream out(outfn.c_str());
	CGAL::output_surface_facets_to_off (out, c2t3);
	return 0;
}

