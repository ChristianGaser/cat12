/*
 * partial_smooth_roi.c  —  Laplacian smoothing limited to a vertex ROI
 *
 *   Updated July 2025 so that it no longer crashes when the input mesh
 *   stores vertices as *single* or faces as *int32/uint32* arrays.
 *   Accepted classes now:
 *       • vertices : double or single, size [N×3]
 *       • faces    : double | single | int32 | uint32, size [M×3]
 *
 *   If an unsupported class is supplied, the MEX throws an informative
 *   error instead of seg‑faulting.
 *
 *   Usage
 *   -----
 *     S2 = partial_smooth_roi(S, ROI [, iterations] [, lambda]);
 *
 *       S.vertices   – [N×3] double|single
 *       S.faces      – [M×3] double|single|int32|uint32 (1‑based)
 *       ROI          – logical|numeric [N×1] mask of vertices to move
 *       iterations   – integer ≥1 (default 10)
 *       lambda       – relaxation 0<λ≤1 (default 0.5)
 *
 *   Compilation (tested R2024a macOS / Apple‑silicon):
 *       mex -largeArrayDims partial_smooth_roi.c
 *
 *   Author: ChatGPT‑o3  |  MIT licence
 */
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef MAX
# define MAX(a,b) ((a)>(b)?(a):(b))
#endif

/* ===================  ADJACENCY LIST  =================================== */
typedef struct {
    mwSize  *idx;   /* concatenated neighbour indices            */
    mwSize  *head;  /* head[i] = start offset of vertex‑i list   */
    mwSize  *deg;   /* deg[i]  = number of neighbours            */
} adjacency_t;

static void free_adjacency(adjacency_t *A)
{
    if (A->idx)  mxFree(A->idx);
    if (A->head) mxFree(A->head);
    if (A->deg)  mxFree(A->deg);
    memset(A,0,sizeof(*A));
}

/* Utility: fetch 0‑based index from faces matrix at linear offset k */
static mwSize get_face_idx(const void *F, mxClassID c, mwSize k)
{
    switch (c) {
        case mxDOUBLE_CLASS:  return (mwSize)(((double  *)F)[k]) - 1U;
        case mxSINGLE_CLASS:  return (mwSize)(((float   *)F)[k]) - 1U;
        case mxINT32_CLASS:   return (mwSize)(((int32_T *)F)[k]) - 1U;
        case mxUINT32_CLASS:  return (mwSize)(((uint32_T*)F)[k]) - 1U;
        default:              return (mwSize)-1; /* should never happen */
    }
}

static void build_adjacency(const mxArray *facesPr, mwSize nV, adjacency_t *A)
{
    const mwSize nF      = mxGetM(facesPr);        /* faces are [M×3] */
    const void  *F       = mxGetData(facesPr);
    const mxClassID cls  = mxGetClassID(facesPr);

    /* ---------- allocate degree & head arrays ---------- */
    A->deg  = (mwSize*)mxCalloc(nV,     sizeof(mwSize));
    A->head = (mwSize*)mxCalloc(nV + 1, sizeof(mwSize));

    /* ---------- first pass: degree count --------------- */
    for (mwSize f = 0; f < nF; ++f) {
        mwSize off = f;   /* row index */
        mwSize v0 = get_face_idx(F, cls, off);
        mwSize v1 = get_face_idx(F, cls, off + nF);
        mwSize v2 = get_face_idx(F, cls, off + 2*nF);
        if (v0>=nV || v1>=nV || v2>=nV)
            mexErrMsgIdAndTxt("partial_smooth_roi:badFace",
                              "Face index out of bounds (are faces 1‑based?)");
        A->deg[v0]+=2; A->deg[v1]+=2; A->deg[v2]+=2;
    }

    /* ---------- prefix sum for heads ------------------- */
    mwSize total = 0;
    for (mwSize i=0;i<nV;++i) {
        A->head[i] = total;
        total += A->deg[i];
    }
    A->head[nV] = total;
    A->idx = (mwSize*)mxCalloc(total, sizeof(mwSize));

    /* reuse deg[] as write‑cursor */
    memset(A->deg, 0, nV*sizeof(mwSize));

    /* ---------- second pass: fill neighbour indices ---- */
    for (mwSize f = 0; f < nF; ++f) {
        mwSize off = f;
        mwSize v0 = get_face_idx(F, cls, off);
        mwSize v1 = get_face_idx(F, cls, off + nF);
        mwSize v2 = get_face_idx(F, cls, off + 2*nF);

        /* v0 <‑‑> v1 */
        A->idx[A->head[v0] + A->deg[v0]++] = v1;
        A->idx[A->head[v1] + A->deg[v1]++] = v0;
        /* v1 <‑‑> v2 */
        A->idx[A->head[v1] + A->deg[v1]++] = v2;
        A->idx[A->head[v2] + A->deg[v2]++] = v1;
        /* v2 <‑‑> v0 */
        A->idx[A->head[v2] + A->deg[v2]++] = v0;
        A->idx[A->head[v0] + A->deg[v0]++] = v2;
    }
}

/* ===================  GATEWAY =========================================== */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2)
        mexErrMsgIdAndTxt("partial_smooth_roi:nargin",
                          "Usage: S2 = partial_smooth_roi(S, ROI [,iter] [,lambda])");

    /* ---- Input structure S ---- */
    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("partial_smooth_roi:invalidS","First input must be a struct.");
    const mxArray *vField = mxGetField(prhs[0],0,"vertices");
    const mxArray *fField = mxGetField(prhs[0],0,"faces");
    if (!vField || !fField)
        mexErrMsgIdAndTxt("partial_smooth_roi:missingFields",
                          "Struct S needs fields 'vertices' and 'faces'.");

    /* ---- validate vertex class ---- */
    const mxClassID vCls = mxGetClassID(vField);
    if (vCls!=mxDOUBLE_CLASS && vCls!=mxSINGLE_CLASS)
        mexErrMsgIdAndTxt("partial_smooth_roi:vertClass",
          "Vertices must be double or single.");

    /* ---- validate face class ---- */
    const mxClassID fCls = mxGetClassID(fField);
    if (fCls!=mxDOUBLE_CLASS && fCls!=mxSINGLE_CLASS &&
        fCls!=mxINT32_CLASS  && fCls!=mxUINT32_CLASS)
        mexErrMsgIdAndTxt("partial_smooth_roi:faceClass",
          "Faces must be double, single, int32 or uint32.");

    /* ---- geometry sizes ---- */
    const mwSize nV = mxGetM(vField);
    if (mxGetN(vField)!=3)
        mexErrMsgIdAndTxt("partial_smooth_roi:vertDim","Vertices must be [N×3].");

    if (mxGetM(prhs[1])!=nV)
        mexErrMsgIdAndTxt("partial_smooth_roi:roiSize",
                          "ROI must have %llu rows (one per vertex)", (unsigned long long)nV);

    /* ---- parameters ---- */
    int    iterations = (nrhs>=3) ? (int)mxGetScalar(prhs[2]) : 10;
    if (iterations < 1) iterations = 1;
    double lambda     = (nrhs>=4) ? mxGetScalar(prhs[3]) : 0.5;
    if (lambda<=0 || lambda>1) lambda = 0.5;

    /* ---- build adjacency once ---- */
    adjacency_t A = {0};
    build_adjacency(fField, nV, &A);

    /* ---- working buffers (double) ---- */
    double *curr = (double*)mxMalloc(nV*3*sizeof(double));
    double *next = (double*)mxMalloc(nV*3*sizeof(double));

    /* copy vertices into curr (promote to double if necessary) */
    if (vCls==mxDOUBLE_CLASS) {
        memcpy(curr, mxGetPr(vField), nV*3*sizeof(double));
    } else { /* single -> double */
        const float *Vs = (float*)mxGetData(vField);
        for (mwSize i=0;i<nV*3;++i) curr[i] = (double)Vs[i];
    }

    /* ---- ROI logical vector ---- */
    mxLogical *roiLog;
    if (mxIsLogical(prhs[1])) {
        roiLog = mxGetLogicals(prhs[1]);
    } else {
        roiLog = (mxLogical*)mxCalloc(nV,sizeof(mxLogical));
        const void *R = mxGetData(prhs[1]);
        const mxClassID rCls = mxGetClassID(prhs[1]);
        for (mwSize i=0;i<nV;++i) {
            double val;
            switch(rCls){
                case mxDOUBLE_CLASS: val=((double*)R)[i]; break;
                case mxSINGLE_CLASS: val=((float *)R)[i]; break;
                case mxINT32_CLASS:  val=((int32_T*)R)[i]; break;
                case mxUINT32_CLASS: val=((uint32_T*)R)[i]; break;
                default:             val=0.0;            break;
            }
            roiLog[i] = (val!=0.0);
        }
    }

    /* ===================  LAPLACIAN SMOOTH  ============================ */
    for (int it=0; it<iterations; ++it) {
        #ifdef _OPENMP
        #pragma omp parallel for if(nV>10000)
        #endif
        for (mwSize v=0; v<nV; ++v) {
            if (!roiLog[v]) { /* untouched vertex */
                next[v]       = curr[v];
                next[v+nV]    = curr[v+nV];
                next[v+2*nV]  = curr[v+2*nV];
                continue;
            }
            mwSize start = A.head[v];
            mwSize deg   = A.deg[v];
            if (deg==0) { /* isolated */
                next[v]       = curr[v];
                next[v+nV]    = curr[v+nV];
                next[v+2*nV]  = curr[v+2*nV];
                continue;
            }
            double sx=0, sy=0, sz=0;
            for (mwSize k=0;k<deg;++k) {
                mwSize nb = A.idx[start+k];
                sx += curr[nb];
                sy += curr[nb+nV];
                sz += curr[nb+2*nV];
            }
            const double invd = 1.0/((double)deg);
            sx *= invd; sy *= invd; sz *= invd;

            next[v]       = curr[v]      + lambda*(sx - curr[v]);
            next[v+nV]    = curr[v+nV]   + lambda*(sy - curr[v+nV]);
            next[v+2*nV]  = curr[v+2*nV] + lambda*(sz - curr[v+2*nV]);
        }
        double *tmp = curr; curr = next; next = tmp; /* swap */
    }

    /* ===================  BUILD OUTPUT ================================ */
    plhs[0] = mxDuplicateArray(prhs[0]); /* shallow copy of input struct */

    mxArray *vOut;
    if (vCls==mxDOUBLE_CLASS) {
        vOut = mxCreateDoubleMatrix(nV,3,mxREAL);
        memcpy(mxGetPr(vOut), curr, nV*3*sizeof(double));
    } else { /* single output */
        vOut = mxCreateNumericMatrix(nV,3,mxSINGLE_CLASS,mxREAL);
        float *Vsingle = (float*)mxGetData(vOut);
        for (mwSize i=0;i<nV*3;++i) Vsingle[i] = (float)curr[i];
    }
    mxSetField(plhs[0],0,"vertices",vOut);

    /* ===================  CLEANUP  ==================================== */
    free_adjacency(&A);
    mxFree(curr); mxFree(next);
    if (!mxIsLogical(prhs[1])) mxFree(roiLog);
}
/* End of file */
