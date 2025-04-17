# Makefile to build make_tetmesh_refine_sphere program
# must first issue 'module load cgal' command
# --- macros
CXX=icc
CXXFLAGS= -DCGAL_DO_NOT_USE_MPZF $(CGAL_INCLUDE)

OBJECTS= make_tetmesh_refine_sphere.o vector.o
LIBS= $(CGAL_LIB)

# --- targets
all:    make_tetmesh_refine_sphere
make_tetmesh_refine_sphere:   $(OBJECTS) 
	$(CXX) $(CXXFLAGS) -o make_tetmesh_refine_sphere   $(OBJECTS) $(LIBS)
        
.SUFFIXES: .cpp 
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<	


# --- remove binary and executable files
clean:
	rm -f make_tetmesh_refine_sphere  $(OBJECTS)


