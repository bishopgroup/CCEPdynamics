/* NORMALIZECOLS.C Normalize the columns of a matrix
 * Syntax: B = normalizecols(A)
 * or B = normalizecols(A,p)
 * The columns of matrix A are normalized so that norm(B(:,n),p) = 1. */
#include <math.h>
#include "mex.h"

/*
% Compute particular solution
PHIP = zeros(param.N1, param.N2, param.N3);
for n = 1:param.N3
    PHIP(:,:,n) = sum(bsxfun(@times,Dn(:,:,:,n),RHOx3(:,:,:)),3);
end

 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Macros for the ouput and input arguments */
    #define P_OUT plhs[0]
    #define D_IN prhs[0]
    #define R_IN prhs[1]
    #define I3_IN prhs[2]
    double *D, *Rr, *Ri, *Pr, *Pi, *I3;
    int iD,iR,iP,n0,n1,n2,n3, NI3;
    size_t KD,KR;
    const mwSize *ND, *NR;
    
    if(nrhs != 3) { /* Check the number of arguments */
        mexErrMsgTxt("Wrong number of input arguments.");
        return;
    }
    
    KD = mxGetNumberOfDimensions(D_IN);
    ND = mxGetDimensions(D_IN);
    D = mxGetPr(D_IN); /* Get the pointer to the data of D */

    KR = mxGetNumberOfDimensions(R_IN);
    NR = mxGetDimensions(R_IN);
    Rr = mxGetPr(R_IN); /* Get the pointer to the data of R */
    Ri = mxGetPi(R_IN); /* Get the pointer to the data of R */
    
    NI3 = mxGetN(I3_IN);
    I3 = mxGetPr(I3_IN);

    if (KD != 4 || KR != 3) {
        mexErrMsgTxt("Wrong size of input arguments.");
        return;
    }
    if (ND[0] != NR[0] || ND[1] != NR[1] || ND[2] != NR[2] || ND[3] != ND[2]) {
        mexErrMsgTxt("Wrong size of input arguments.");
        return;
    }

    P_OUT = mxCreateNumericArray(KR, NR, mxDOUBLE_CLASS, mxCOMPLEX);
    Pr = mxGetPr(P_OUT); 
    Pi = mxGetPi(P_OUT); 
    
    for(n0 = 0; n0 < ND[0]; n0++) {
        for(n1 = 0; n1 < ND[1]; n1++) {
            for(n2 = 0; n2 < ND[2]; n2++) {
                iP = n0 + ND[0]*(n1 + ND[1]*n2);
                Pr[iP] = 0;
                Pi[iP] = 0;
                for(n3 = 0; n3 < NI3; n3++) {
                    iR = n0 + ND[0]*(n1 + ND[1]*(I3[n3]-1));
                    iD = n0 + ND[0]*(n1 + ND[1]*((I3[n3]-1) + ND[2]*n2));
                    Pr[iP] += D[iD]*Rr[iR];
                    Pi[iP] += D[iD]*Ri[iR];
                }
            }
        }
    }
    
    return;
}