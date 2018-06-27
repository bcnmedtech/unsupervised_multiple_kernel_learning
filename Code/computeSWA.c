/*========================================================*
 *  This Unsupervised Multiple Kernel Learning (U-MKL) package is (c) BCNMedTech
 *  UNIVERSITAT POMPEU FABRA
 * 
 *  This U-MKL package is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 * 
 *  This U-MKL package is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 * 
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 *  Author:
 *  Sergio Sanchez-Martinez 
 *
 *  Contributors: 
 *  Nicolas Duchateau
 *  Gemma Piella
 *  Constantine Butakoff
 *========================================================*/

/*========================================================*
 * computeSWA.c
 *
 * The calling syntax from Matlab is:
 *
 *		[S_W_A , S_D_A] = computeSWA(K_tot_ALL, A, W, Diag);
 *
 * This is a MEX file for MATLAB.
 *========================================================*/

#include "mex.h"

/* The computational routine */
void computeSWA(double *K_tot_ALL, double *A, double *W, double *Diag, double *S_W_A, double *S_D_A, int N, int NA, int m, int i1, int i2)
{
    int r,l,i,j,c,c1,c2,k;

    double *tmp;
    tmp = (double*)mxMalloc(sizeof(double)*m*NA);
    
    double *KiA;
    KiA = (double*)mxMalloc(sizeof(double)*N*m*NA);

    double tmpV;
    
    for(i=0; i<N; i++)
    {
        for(c=0; c<m; c++)
        {
            for(r=0; r<NA; r++)
            {
                tmpV = 0.;
                for(l=0; l<N; l++)
                {
                    tmpV += K_tot_ALL[i + l*N + c*N*N] * A[l + r*N];
                }
                KiA[i + N*c + r*N*m] = tmpV;
            }
        }
    }
    
    for(i=i1; i<i2; i++)
    {
        for(j=(i+1); j<N; j++)
        {
            double Wij = W[i + j*N];
            
            if(Wij != 0.)
            {
                for(c1=0; c1<m; c1++)
                {
                    for(c2=c1; c2<m; c2++) 
                    {
                        tmpV = 0.;
                        for(k=0; k<NA; k++)
                        {
                            double tmp1 = KiA[i + N*c1 + k*N*m] - KiA[j + N*c1 + k*N*m];
                            double tmp2 = KiA[i + N*c2 + k*N*m] - KiA[j + N*c2 + k*N*m];
                            tmpV += tmp1 * tmp2;
                        }
                        S_W_A[c1 + c2*m] += 2. * Wij * tmpV;
                    }
                }
            }
        }
                
        
        for(c=0; c<m; c++)
        {
            for(r=0; r<NA; r++)
            {
                tmpV = 0.;
                for(l=0; l<N; l++)
                {
                    tmpV += K_tot_ALL[i + l*N + c*N*N] * A[l + r*N];
                }
                tmp[c + r*m] = tmpV;
            }
        }        
        
        const double Diagii = Diag[i];

        for(c1=0; c1<m; c1++)
        {
            for(c2=c1; c2<m; c2++) 
            {
                tmpV = 0.;
                for(k=0; k<NA; k++)
                {
                    tmpV += (tmp[c1 + k*m] * tmp[c2 + k*m]);
                }
                S_D_A[c1 + c2*m] += Diagii * tmpV;
            }
        }        
    }
    mxFree(tmp);
    mxFree(KiA);
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 6)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "6 inputs required.");
    }
    if(nlhs != 2)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "2 outputs required.");
    }
    
    int nDimNum = mxGetNumberOfDimensions(prhs[0]);
    const int *pDims = mxGetDimensions(prhs[0]);  

    int nDimNumA = mxGetNumberOfDimensions(prhs[1]);
    const int *pDimsA = mxGetDimensions(prhs[1]);  
    int N = pDims[0];
    int NA = pDimsA[1];
    int m = pDims[2];
    
    double *K_tot_ALL = mxGetPr(prhs[0]); 
    double *A = mxGetPr(prhs[1]); 
    double *W = mxGetPr(prhs[2]); 
    double *Diag = mxGetPr(prhs[3]); 

    double *i1 = mxGetPr(prhs[4]);
    double *i2 = mxGetPr(prhs[5]);

    
    /*mexPrintf("ComputeSWB:Requested range: %d-%d\n",(int)i1[0],(int)i2[0]);*/
    
    
    int outDims[2];
    outDims[0] = m;
    outDims[1] = m;
    
    plhs[0] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    double *S_W_A = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    double *S_D_A = mxGetPr(plhs[1]);

    computeSWA(K_tot_ALL,A,W,Diag,S_W_A,S_D_A,N,NA,m,(int)i1[0],(int)i2[0]);
}

