// Matlab usage:
// r = intersect_seg_triangles_mex(p1,p2,te,p,tiny);
// p1, p2 define the segment (1x3 vectors of double)
// te, p define the triangles:
//   te: nx3
//   p : mx3
// tiny: is the epsilon value used for round-off errors in floating point calculations
// 
// returns non-zero r if a segment is intersecting one of the triangles in 'te'


/* How to mex
*	if your compiler supports OpenMP, use this:
*   GCC:
*   mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" intersect_seg_triangles_mex.cpp geomath.cpp vector.cpp
*
*	Intel/{Mac,Linux}
*   mex -v CXXFLAGS="\$CXXFLAGS -openmp" LDFLAGS="\$LDFLAGS -openmp" intersect_seg_triangles_mex.cpp geomath.cpp vector.cpp
*
*   Microsoft Visual C++ 2008
*   mex -v COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp" intersect_seg_triangles_mex.cpp geomath.cpp vector.cpp
*
*	Intel/{Windows}
*	mex -v CXXFLAGS="\$CXXFLAGS -Qopenmp" LDFLAGS="\$LDFLAGS -Qopenmp" intersect_seg_triangles_mex.cpp geomath.cpp vector.cpp
*
*	If it doesn't support, then just use:
*	mex -v intersect_seg_triangles_mex.cpp geomath.cpp vector.cpp
*/

#include "geomath.h"
#include "mex.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define points(i,j) points[(i)+np*(j)]
#define ele(i,j) ele[(i)+ne*(j)]
#define t(i,j) t[(i)+ne*(j)]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs!=1 || nrhs!=5)
        mexErrMsgTxt("intersect_seg_triangles_mex: needs 5 input and 1 output\n");
    
    unsigned long np = mxGetM(prhs[3]);
    unsigned long ne = mxGetM(prhs[2]);
    
    //mexPrintf(" ne: %d, np: %d\n", ne, np);
    
    plhs[0] = mxCreateNumericMatrix(1, 1, mxINT8_CLASS, mxREAL);
    
    double *p1 = mxGetPr(prhs[0]);
    double *p2 = mxGetPr(prhs[1]);
    double eps = *(mxGetPr(prhs[4]));
    //mexPrintf(" eps = %f\n", eps);
    
    double *t;
    ULONG *ele = new ULONG[ne*3];
    if (mxIsDouble(prhs[2])) {
        t = mxGetPr(prhs[2]);
        for (ULONG i=0; i<ne; ++i) {
            for (int j=0; j<3; ++j)
            {
                ele(i,j) = (ULONG) t(i,j);
                //mexPrintf("ele(i,j) = %d\n", ele(i,j));
            }
        }
    }
    else
        mexErrMsgTxt("intersect_seg_triangles_mex: input element list should be in double\n");

    //mexPrintf(" done double to ulong\n");
    
	double *points = mxGetPr(prhs[3]);
	
    int st = 0;
    bool flag = false;

    
    #ifdef _OPENMP
	omp_set_num_threads(omp_get_num_procs());
// 	#pragma omp single
// 	mexPrintf(" CPUs available: %d (%d)\n",omp_get_num_procs(), omp_get_thread_num());
    #endif
    
    int i,j,ii;
    double I[3], tp1[3], tp2[3], tp3[3];
    #pragma omp parallel default(none) \
        private(i,j,ii,tp1,tp2,tp3,I) \
        firstprivate(st) \
        shared(ne,flag,points,ele,eps,p1,p2,np,std::cout)
    {
        #pragma omp for
        for (i=0; i<ne; ++i)
        {
            for (j=0; j<3; ++j)
            {
                //mexPrintf(" component j: %d\n", j);
                tp1[j] = points(ele(i,0)-1, j);
                tp2[j] = points(ele(i,1)-1, j);
                tp3[j] = points(ele(i,2)-1, j);
            }
            //mexPrintf("calling intersect: %d\n",i);
            st = intersect_RayTriangle(p1, p2, tp1, tp2, tp3, I, eps);
            //mexPrintf(" st is %d\n", st);
            if (st == 1 || st > 300 || (st<=12 && st>=10))
            {
                flag = true;
                mexPrintf(" INTERSECTIOIN: (on %d)\n", st);
                #pragma omp strict
                {
                    std::cout << " st is: " << st << std::endl;
                    std::cout.flush();
                }
                mexPrintf("p1=[%.12f, %.12f, %.12f];\n", p1[0],p1[1],p1[2]);
                mexPrintf("p2=[%.12f, %.12f, %.12f];\n", p2[0],p2[1],p2[2]);
                mexPrintf(" element: %d\n\n", i+1);
                for (ii=0; ii<3; ++ii)
                {
                    mexPrintf("%.12f %.12f %.12f\n",points(ele(i,ii)-1, 0),points(ele(i,ii)-1, 1),points(ele(i,ii)-1, 2));
                }
            }
            else if (st<=22 && st>=20)
            {
                // mexPrintf(" got 20!!! values\n");
                // TODO: need to finish this part
            }
            else if (st<=202 && st>=200)
            {
                // mexPrintf(" got 200!!! values\n");
            }
        }
    }
    
    char *c = (char *) mxGetData(plhs[0]);
    if (flag)
        c[0] = 1;
    else
        c[0] = 0;
    
    delete [] ele;
}