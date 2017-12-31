#include "mex.h"
#include <math.h>
#include "dg_chamfer.h"

/*

This is a modified version of mm_geodesic_voronoi.c where you can add
geometric constraints during the front propagation :
 - increase distance to "jump" a sulcus
 - impose a seed to be at a certain place

*/
#define _DG_DEBUG 0
#undef _DG_DEBUG1
#define MAGIC_NUMBER 1

/* --- GEODESIC VORONOI SUBROUTINE --- */
void geodesic_voronoi(double *img,         /* domain image mask */
			   		  unsigned int *size,  /* size of domain image */
			   		  double *voxsize,     /* size of voxels (in mm) */
					  double *seeds,       /* positions of the seeds */
					  int nbseeds,         /* number of seeds*/
					  double *aniso,       /* Anisotropic ponderation */
					  double *vor,         /* voronoi diagram */
					  double *dmap,        /* geodesic distance map */
					  char *dist) {        /* type of distance */
	int i, j, a;
	int k = 1, l = 1, p = 1, s = 1;
	unsigned int m = size[0], n = size[1], o = size[2];
	
	mask *cmask = NULL;
	queue q;
	int notempty = 1;
	int current_buck = 0;
	int x, y, z;
	double newd;
	
	/*mexPrintf("Ca y est je rentre ...");*/
	
	/* --- CHAMFER MASK INITIALIZATION --- */
	if (!strcmp(dist,"d4")) {
		cmask = maskd4;
		k = l = p = 1;
		s = 6;
	}
	else if (!strcmp(dist,"d8")) {
		cmask = maskd8;
		k = l = p = 1;
		s = 18;
	}
	else if (!strcmp(dist,"d34")) {
		cmask = maskd34;
		k = l = p = 1;
		s = 26;
	}
	else if (!strcmp(dist,"d5711")) {
		cmask = maskd5711;
		k = l = p = 2;
		s = 100;
		mexWarnMsgTxt("Geodesic distances may be underestimated with d5711.");
	}
	else {
		mexErrMsgTxt("Unknown distance type.");
	}
	
	/* --- VOXEL SIZE WEIGHTING --- */
	
	
	/* --- BUCKET, DISTMAP AND VORONOI INITIALIZATION --- */
	q.b = mxMalloc(sizeof(bucket));
	q.nbbucket = 1;
	
	/* distmap and voronoi imitialization */
	for (i=0;i<m*n*o;i++)
		if (img[i] > 0.0) {
			dmap[i] = mxGetInf();  /* inside domain */
			vor[i]  = mxGetNaN();
		}
		else {
			dmap[i] = -1; /* outside domain */
			vor[i]  = -1;
		}
	
	mexPrintf("Avant Init Nb seed : %d\n", nbseeds);	
		
	/* bucket initialization and distmap and voronoi update */
	q.b[0].size = nbseeds;
	q.b[0].pt = mxMalloc(nbseeds*sizeof(point3d));
	q.b[0].memsize = nbseeds;
	for (i=0;i<nbseeds;i++) {
		q.b[0].pt[i].x = (int)(seeds[i]+0.5) - 1;
		q.b[0].pt[i].y = (int)(seeds[i + nbseeds]+0.5) - 1;
		q.b[0].pt[i].z = (int)(seeds[i + 2 * nbseeds]+0.5) - 1;
		
		/*mexPrintf("Indice %d : %d %d %d \n",i, q.b[0].pt[i].x , q.b[0].pt[i].y , q.b[0].pt[i].z); */
		
		dmap[q.b[0].pt[i].x+m*q.b[0].pt[i].y+m*n*q.b[0].pt[i].z] = 0.0;
		vor[q.b[0].pt[i].x+m*q.b[0].pt[i].y+m*n*q.b[0].pt[i].z]  = i+1;
	}
	
	/*mexPrintf("Apres Init \n!");*/
	
	/* --- BUCKET SORTING --- */
	while (notempty) {
#if _DG_DEBUG
		mexPrintf("\r%d",current_buck);fflush(NULL);
#endif
		/* loop over the not empty bucket with the smallest value */
		for (i=0;i<q.b[current_buck].size;i++) {
			/* loop over the chamfer mask */
			for (a=0;a<s;a++) {
				x = q.b[current_buck].pt[i].x + cmask[a].x;
				y = q.b[current_buck].pt[i].y + cmask[a].y;
				z = q.b[current_buck].pt[i].z + cmask[a].z;
				/* point inside image ? */
				if ((x >= 0) && (y >= 0) && (z >= 0)
					 && (x < m) && (y < n) && (z < o)) {
#ifdef _DG_DEBUG1
				  /*mexPrintf("%d : %d %d %d ; %f\n",a, x, y, z,dmap[x+m*y+m*n*z]);*/
#endif
					/* point inside domain ? */
					if (dmap[x+m*y+m*n*z] > 0.0) {
#ifdef _DG_DEBUG
						mexPrintf("%d %d %d  ",x,y,z);
#endif
						/* Weighted distances */
						if (aniso) {
							newd = dmap[q.b[current_buck].pt[i].x
								  	 + m*q.b[current_buck].pt[i].y
								   	 + m*n*q.b[current_buck].pt[i].z]
							   		 + aniso[x+m*y+m*n*z] * cmask[a].v;
						
							/* connexite conditionnelle */
							/* 6 connexite (deplacement de 1)*/
							/*if (abs(cmask[a].x) + abs(cmask[a].y) + abs(cmask[a].z) == 1) {
								newd = dmap[q.b[current_buck].pt[i].x
								  	 + m*q.b[current_buck].pt[i].y
								   	 + m*n*q.b[current_buck].pt[i].z]
							   		 + aniso[x+m*y+m*n*z] * cmask[a].v;
							}*/
							/* 18 connexite (deplacement de 2) */
							/*else if (abs(cmask[a].x) + abs(cmask[a].y) + abs(cmask[a].z) == 2){
								newd = dmap[q.b[current_buck].pt[i].x
								  	 + m*q.b[current_buck].pt[i].y
								   	 + m*n*q.b[current_buck].pt[i].z]
							   		 + aniso[x+m*y+m*n*z] * cmask[a].v;
							}*/
							/* 26 connexite (deplacement de 3) */
							/*else {
								newd = dmap[q.b[current_buck].pt[i].x
								  	 + m*q.b[current_buck].pt[i].y
								   	 + m*n*q.b[current_buck].pt[i].z]
							   		 + aniso[x+m*y+m*n*z] * cmask[a].v;
							}*/
						}
						/* Non-weighted distances */
						else
							newd = dmap[q.b[current_buck].pt[i].x
									+ m*q.b[current_buck].pt[i].y
									+ m*n*q.b[current_buck].pt[i].z]
									+ cmask[a].v;
#ifdef _DG_DEBUG
						mexPrintf("%f   %f  ",newd, dmap[x+m*y+m*n*z]);
#endif
						/* shorter distance ? */
						if (newd < dmap[x+m*y+m*n*z]) {
#ifdef _DG_DEBUG
							mexPrintf("   add");
#endif
							/* new value for that distance */
							dmap[x+m*y+m*n*z] = newd;
							/* Update Voronoi diagram */
							vor[x+m*y+m*n*z] = vor[q.b[current_buck].pt[i].x
								    + m*q.b[current_buck].pt[i].y
								    + m*n*q.b[current_buck].pt[i].z];
						 	/* add in bucket corresponding to newd */
							/* queue size update */
							if (q.nbbucket <= (int)newd) {
								q.b = mxRealloc(q.b,((int)newd+1) * sizeof(bucket));
								for (j=q.nbbucket;j<(int)newd+1;j++) {
									q.b[j].size = 0;
									q.b[j].memsize = 2;
									q.b[j].pt = mxMalloc(2*sizeof(point3d));
								}
								q.nbbucket = (int)newd+1;
							}
							/* bucket size update */
							if (q.b[(int)newd].size >= q.b[(int)newd].memsize) {
								q.b[(int)newd].pt = mxRealloc(q.b[(int)newd].pt,q.b[(int)newd].memsize*2*sizeof(point3d));
								q.b[(int)newd].memsize *=2;
							}	
							q.b[(int)newd].pt[q.b[(int)newd].size].x = x;
							q.b[(int)newd].pt[q.b[(int)newd].size].y = y;
							q.b[(int)newd].pt[q.b[(int)newd].size].z = z;
							q.b[(int)newd].size++;
					  	}
					}
				}
			}
		}
#ifdef _DG_DEBUG
		/*mexPrintf("\nNombre d'objets : %d\n",q.b[current_buck].size);*/
		/*mexPrintf("Allocation      : %d\n",q.b[current_buck].memsize);*/
#endif
		mxFree(q.b[current_buck].pt);
		/* empty buckets ? */
		if (q.nbbucket == current_buck + 1)
			notempty = 0;
		else
			current_buck++;
	}
#ifdef _DG_DEBUG
	mexPrintf("\b                      \r");
#endif
			
	return;
}

/* --- GATEWAY FUNCTION --- */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) {

double *img = NULL, *aniso = NULL, *distmap = NULL, *vor = NULL;
double *voxsize = NULL, *seeds = NULL;
unsigned int size[3], sdist, nbseeds;
char *dist = NULL;
mxArray *dmap = NULL, *vord = NULL;

/* Check for proper number of arguments. */
if (nrhs == 0)
	mexErrMsgTxt("Not enough input arguments.");
else if (nrhs > 5)
	mexErrMsgTxt("Too many input arguments.");
else if (nlhs > 2)
	mexErrMsgTxt("Only one output argument allowed.");

/* The input IMG must be a real double array. */
/*mexPrintf("dim : %d\n",mxGetNumberOfDimensions(prhs[0]));*/
if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
	!(mxGetNumberOfDimensions(prhs[0])==2 || 
	mxGetNumberOfDimensions(prhs[0])==3))
	mexErrMsgTxt("Input img must be a real double array.");
size[0] = mxGetDimensions(prhs[0])[0];
size[1] = mxGetDimensions(prhs[0])[1];
if (mxGetNumberOfDimensions(prhs[0])==3)
	size[2] = mxGetDimensions(prhs[0])[2];
else 
	size[2] = 1;
img = mxGetPr(prhs[0]);

/* The input VOXSIZE must be a [x y z] array */
if (nrhs > 1) {
	if (!mxIsDouble(prhs[1]) 
			|| mxIsComplex(prhs[1])
			|| mxGetNumberOfElements(prhs[1])!=3)
	mexErrMsgTxt("Input voxsize must be a 3 elements vector.");
	voxsize = mxGetPr(prhs[1]);
}
else {
	voxsize = mxMalloc(3*sizeof(double));
	voxsize[0] = voxsize[1] = voxsize[2] = 1.0;
}

/* The input SEEDS must be a real double array [n x 3] */
if (!mxIsDouble(prhs[2]) || 
	 mxIsComplex(prhs[2]) || 
	 mxGetNumberOfDimensions(prhs[2]) != 2 ||
	 mxGetDimensions(prhs[2])[1] != 3)
	mexErrMsgTxt("Input seeds must be a real double array [n x 3].");
seeds = mxGetPr(prhs[2]);
nbseeds = mxGetDimensions(prhs[2])[0];
/*mexPrintf("number of seeds : %d\n",nbseeds);*/


/* The input DIST must be a string */
if (nrhs > 3) {
	if (!mxIsChar(prhs[3]))
		mexErrMsgTxt("Input dist must be a string.");
	sdist = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
	dist = mxCalloc(sdist, sizeof(char));
	mxGetString(prhs[3], dist, sdist);
}
else {
	dist = "d34";
}


/* The input ANISO must be a real double array */
if (nrhs > 4) {
	/* mexPrintf("dim : %d\n",mxGetNumberOfDimensions(prhs[4])); */
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
		!(mxGetNumberOfDimensions(prhs[4])==2 || 
		mxGetNumberOfDimensions(prhs[4])==3))
		mexErrMsgTxt("Input aniso must be a real double array.");
	if (mxGetNumberOfDimensions(prhs[0]) != mxGetNumberOfDimensions(prhs[4]))
		mexErrMsgTxt("Inputs img and aniso are of different sizes.");
	if ((mxGetDimensions(prhs[4])[0] != size[0]) 
		|| (mxGetDimensions(prhs[4])[1] != size[1]))
		mexErrMsgTxt("Inputs img and aniso are of different sizes.");
	aniso = mxGetPr(prhs[4]);
	/* mexPrintf("Anisotropic front propagation\n"); */
}
else {
	/* mexPrintf("Isotropic front propagation\n"); */
}

/* Create mxArray's */
dmap = mxCreateNumericArray(3, (const int*)size, mxDOUBLE_CLASS, mxREAL);
distmap = mxGetPr(dmap);

vord = mxCreateNumericArray(3, (const int*)size, mxDOUBLE_CLASS, mxREAL);
vor = mxGetPr(vord);

/* Call the Geodesic Voronoi subroutine. */
geodesic_voronoi(img, size, voxsize, seeds, nbseeds, aniso, vor, distmap, dist);

/* Assign pointers to output */
if (nlhs > 0) plhs[0] = vord;
else mxDestroyArray(vord);

if (nlhs > 1) plhs[1] = dmap;
else mxDestroyArray(dmap);

/* Free memory */
if (nrhs < 2) mxFree(voxsize);

}
