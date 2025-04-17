/*********************************************************************
 * Modified 20 oct 2017, Matt McGarry
 * Corrected a number of errors
 * Generalized to interpolate arbitraty tetrahedral nodal values
 * Set the (1,1,1) voxel to be at [x,y,z]=[0,0,0]  
 *
 * Previous version by 
 *  AJ Pattison
 *  Graduate Student
 *  Dartmouth College
 *  May 2008
 *
 *  SHEAR_STACK_MEX.c
 *
 *  Modified from: morph_stack_mex.c
 *
 *  Songbai Ji
 *  Dartmouth College
 *  Image Guided Neurosurgery
 *  April 2008
 ********************************************************************
 C mex function that interpolates nodal tet values back onto a volume.
 tet_stack_mex(stack_size, voxel_size, int32(elm), nod, tetval(:,ii));
	Input: image stack (stack), voxel size (vsize), element matrix (el), xyz, tetvals
	Output: Image stack of tet.
	tetbrain = tet_stack_mex(stack, stack_size, (int32)elem, nod, tet); 
	
    Note: should make sure elem is of 32-bit integer data type, so that the program can run 
		on both 32-bit or 64-bit machines. 
*/

#include "mex.h"
#include <math.h>

/* #define p(i,j,k) pr[(i-1) + (j-1)*M + (k-1)*M*N ]*/  /* this is 1-based, like Matlab */
/* #define p(i,j,k) pr[(i) + (j)*M + (k)*M*N ]  */ /* this is 0-based, like in C */

#define stack(i,j,k) stack_pr[(i) + (j)*stack_M + (k)*stack_M*stack_N ]
#define stack_tet(i,j,k) stack_tet_pr[(i) + (j)*stack_M + (k)*stack_M*stack_N ]
#define el(i,j) (el_pr[(i) + (j)*nelm ] -1) /* IMPORTANT: mesh elements span from 1 to whatever, therefore, we need 
												   to subtract by 1 to make it 0-based */
#define nd(i,j) nd_pr[(i) + (j)*nnod ]
#define tet(i,j) tet_pr[(i) + (j)*nnod ]
#define nd2(i,j) nd_pr2[(i) + (j)*nnod ]
#define tet2(i,j) tet_pr2[(i) + (j)*nnod ]

#define rstack(i,j,k) rstack_pr[(i) + (j)*stack_M + (k)*stack_M*stack_N ]

bool withinTet(int *el_pr, mwSize nelm, double *nd_pr, mwSize nnod, int ind, double p[3], double basis[4]);
double vol_tet(double a[3], double b[3], double c[3], double d[3]) ;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double *stack_pr, *nd_pr, *tet_pr, *stack_tet; /* stack_disp is the tet stack */
	int *el_pr; /* because this would enter as subscripts for nodal positions, has to be integer */
	mwSize stack_M, stack_N, stack_Z, NUM_VOXELS, *dims;
	mwSize nelm, nnod; /* number of elements and number of nodes */
	double vx, vy, vz; /* voxel sizes in mm along x, y, z, which correspond to column, row and slice direction */
	double minx, miny, minz, maxx, maxy, maxz; /* bounding box for a tet */
	int rx, ry, rz; /* voxels within the bounding box */
	double voxel[3], basis[4]; /* basis[0] is total volume, 1-4 sub volumes */
	double tet_val; /* tet for voxel */
	int i, j, ind;
    double *rstack_pr; /* returned stack */
	double *nd_pr2, *tet_pr2; /* temporary variables to replace input nd and tet, so they are not modified when return to Matlab */
	bool isInside; 

	/* check data type of input */
	if (nrhs != 5) mexErrMsgTxt("5 inputs are required.");

	if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS )
		mexErrMsgTxt("image stack must be double.");	
	if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS )
		mexErrMsgTxt("image stack voxel size must be double.");
	if (mxGetClassID(prhs[2]) != mxINT32_CLASS )
		mexErrMsgTxt("element matrix must be INT32.");
	if (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS )
		mexErrMsgTxt("nodal positions must be double.");
	if (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS )
		mexErrMsgTxt("tets must be double.");

	/* check input further and get input */
	if (mxGetNumberOfDimensions(prhs[0]) == 3) {
		stack_pr = mxGetPr(prhs[0]);
	
		stack_M = mxGetDimensions(prhs[0])[0]; 
		stack_N = mxGetDimensions(prhs[0])[1];
		stack_Z = mxGetDimensions(prhs[0])[2]; 
	} else {
		mexErrMsgTxt("The image stack must be a 3D matrix.");
	}
	
	if (mxGetNumberOfElements(prhs[1]) ==3) {
		vx = mxGetPr(prhs[1])[0]; vy = mxGetPr(prhs[1])[1]; vz = mxGetPr(prhs[1])[2]; 
	} else
		mexErrMsgTxt("Voxel size must be of vector of 3.");

	if (mxGetN(prhs[2]) ==4) {
		nelm = mxGetM(prhs[2]);
		el_pr = mxGetPr(prhs[2]);
	} else
		mexErrMsgTxt("Element input must be 4 columns.");

	if (mxGetN(prhs[3]) ==3) {
		nnod = mxGetM(prhs[3]);
		nd_pr = mxGetPr(prhs[3]);
	} else
		mexErrMsgTxt("Nodal input must be 3 columns.");
	if (mxGetN(prhs[4]) ==1) {
		tet_pr = mxGetPr(prhs[4]);
	} else
		mexErrMsgTxt("tet input must be 1 columns.");

	/* create temp and return arrays */
	dims = mxGetDimensions(prhs[0]); 
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	rstack_pr = (double *) mxGetData(plhs[0]);

	/* replace nd and tet with these temporary variables, so that nd and tet not modified */
	nd_pr2 = mxMalloc(nnod*3*sizeof(double));
	tet_pr2 = mxMalloc(nnod*sizeof(double));
	
	/* check memory allocation */
	if (rstack_pr == NULL || nd_pr2 == NULL || tet_pr2 == NULL)
		mexErrMsgTxt("Error in memory allocation.");

	/* first scale nodal info (m) to voxel unit; here we copied the input matrices to local variable.
	This is because nd and tet may be changed otherwise */
	for (i=0; i<nnod; i++) {
		nd2(i,0) = nd(i,0)/vy +1.0;  
		nd2(i,1) = nd(i,1)/vx +1.0;
		nd2(i,2) = nd(i,2)/vz +1.0;
        tet2(i,0) = tet(i,0);
	}

	/* then, loop through all elements */
	for (i=0; i<nelm; i++) {
		/* find the bounding box */
		minx=999; miny=999; minz=999;
		maxx=-999; maxy=-999; maxz=-999;
		for (j=0; j<4; j++) {
			if (minx>nd2(el(i,j), 0)) minx=nd2(el(i,j), 0);
			if (miny>nd2(el(i,j), 1)) miny=nd2(el(i,j), 1);
			if (minz>nd2(el(i,j), 2)) minz=nd2(el(i,j), 2);
			if (maxx<nd2(el(i,j), 0)) maxx=nd2(el(i,j), 0);
			if (maxy<nd2(el(i,j), 1)) maxy=nd2(el(i,j), 1);
			if (maxz<nd2(el(i,j), 2)) maxz=nd2(el(i,j), 2);
		}
		
		minx=minx-0.51; miny=miny-0.51; minz=minz-0.51;
		maxx=maxx+0.51; maxy=maxy+0.51; maxz=maxz+0.51; /* relax a bit */
		/* here, we should check if minx, maxx, etc are outside the image stack range. However, if a voxel is outside
		this range, it should be outside the mesh as well if nothing unexpected happened. So maybe we don't need to
		explicitly check.  In addition, this is actually checked later (in "isInside" clause). */
		if (minx<1) minx=1; if (miny<1) miny=1; if (minz<1) minz=1;
		if (maxx>stack_N) maxx=stack_N; if (maxy>stack_M) maxy=stack_M; if (maxz>stack_Z) maxz=stack_Z;

		/* now navigate through all enclosed voxels */
		for (rx=(int)minx; rx<=(int)maxx+1; rx++) {
			for (ry=(int)miny; ry<=(int)maxy+1; ry++) {
				for (rz=(int)minz; rz<=(int)maxz+1; rz++) {
					/* ignore voxels outside tet */
					voxel[0]=rx; voxel[1]=ry; voxel[2]=rz; 
					isInside = withinTet(el_pr, nelm, nd_pr2, nnod, i, voxel, basis) ;

					if (isInside){
						tet_val = tet2(el(i,0),0)*basis[1] + tet2(el(i,1),0)*basis[2] + 
							tet2(el(i,2),0)*basis[3] + tet2(el(i,3),0)*basis[4];
						tet_val = tet_val/basis[0];

      					ind = (ry-1) + (rx-1)*stack_M + (rz-1)*stack_M*stack_N;
                        rstack_pr[ind] = tet_val;

					} /* if isInside */
				} /* for rz */
			} /* for ry */
		} /* for rx */
	} /* for (i=0; i<nelm; i++) */ 
	
	mxFree(nd_pr2); mxFree(tet_pr2);
}

bool withinTet(int *el_pr, mwSize nelm, double *nd_pr, mwSize nnod, int ind, double p[3], double basis[4]) {
	int j; 
	double a[3], b[3], c[3], d[3]; /* four nodes to form a tet */
	double vol1, vol2, vol3, vol4, volT; /* volumes for 5 tets */
	
	/* construct 5 tets */
	for (j=0; j<3; j++) a[j] = nd(el(ind,0),j);
	for (j=0; j<3; j++) b[j] = nd(el(ind,1),j);
	for (j=0; j<3; j++) c[j] = nd(el(ind,2),j);
	for (j=0; j<3; j++) d[j] = nd(el(ind,3),j);
	volT = fabs(vol_tet(a,b,c,d));

	vol1 = fabs(vol_tet(p,b,c,d));
	vol2 = fabs(vol_tet(a,p,c,d));
	vol3 = fabs(vol_tet(a,b,p,d));
	vol4 = fabs(vol_tet(a,b,c,p));
    
	basis[0] = volT; basis[1]=vol1; basis[2]=vol2; basis[3]=vol3; basis[4]=vol4;

	if (fabs( vol1+vol2+vol3+vol4 - volT ) < 1e-10)
		return true;
	else
		return false;
}

double vol_tet(double aa[3], double bb[3], double cc[3], double dd[3])
{
	int i;
	double volume, a[3], b[3], c[3]; 

/* we must make sure that aa, bb, cc and dd values are not changed in this function */
	for (i=0; i<3; i++) {
		a[i] = aa[i]-dd[i];
		b[i] = bb[i]-dd[i];
		c[i] = cc[i]-dd[i];
	}
	volume = a[0]*(b[1]*c[2] - b[2]*c[1]) + 
		a[1]*(b[2]*c[0] - b[0]*c[2]) + 
		a[2]*(b[0]*c[1] - b[1]*c[0]) ;
	return volume; 
	/* The real volume should be further divided by 6.
	However, because we are using the ratio of volumes, not really 
	the volume itself, we are fine. */
}
