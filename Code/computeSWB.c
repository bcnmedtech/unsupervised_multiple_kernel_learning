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
 * computeSWB.c
 *
 * The calling syntax from Matlab is:
 *
 *		[S_W_B , S_D_B] = computeSWB(K_tot_ALL, betas, W, Diag);
 *
 * This is a MEX file for MATLAB.
 *========================================================*/

#include "mex.h"

/* The computational routine */
void computeSWB(double *K_tot_ALL, double *betas, double *W, double *Diag, double *S_W_B, double *S_D_B, int N, int m, int i1, int i2)
{
    int r,l,i,j,c;

    
    double *tmp;
    tmp = (double*)mxMalloc(sizeof(double)*N);
    
    double tmpV;
    
    for(i=i1; i<i2; i++)
    {
        for(j=(i+1); j<N; j++)
        {
            double Wij = W[i + j*N];
            
            if(Wij != 0.)
            {
                for(r=0; r<N; r++)
                {
                    tmpV = 0.;
                    for(c=0; c<m; c++)
                    {
                        tmpV += (K_tot_ALL[i + r*N + c*N*N] - K_tot_ALL[j + r*N + c*N*N]) * betas[c];
                    }
                    tmp[r] = tmpV;
                }

                for(r=0; r<N; r++)
                {
                    for(l=r; l<N; l++)
                    {
                        S_W_B[r + l*N] += 2. * Wij * (tmp[r] * tmp[l]);
                    }
                }
            }
        }
        
        
        for(r=0; r<N; r++)
        {
            tmpV = 0.;
            for(c=0; c<m; c++)
            {
                tmpV += K_tot_ALL[i + r*N + c*N*N] * betas[c];                
            }
            tmp[r] = tmpV;
        }        
        
        const double Diagii = Diag[i];

        tmpV = 0.;
        for(r=0; r<N; r++)
        {
            for(l=r; l<N; l++) 
            {
                S_D_B[r + l*N] += Diagii * (tmp[r] * tmp[l]);
            }
        }
    }
    mxFree(tmp);
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
    int N = pDims[0];
    int m = pDims[2];
    
    double *K_tot_ALL = mxGetPr(prhs[0]); 
    double *betas = mxGetPr(prhs[1]); 
    double *W = mxGetPr(prhs[2]); 
    double *Diag = mxGetPr(prhs[3]); 
    double *i1 = mxGetPr(prhs[4]);
    double *i2 = mxGetPr(prhs[5]);

    
    mexPrintf("ComputeSWB:Requested range: %d-%d\n",(int)i1[0],(int)i2[0]);
    
    
    int outDims[2];
    outDims[0] = N;
    outDims[1] = N;
    
    plhs[0] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    double *S_W_B = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    double *S_D_B = mxGetPr(plhs[1]);
        
    computeSWB(K_tot_ALL,betas,W,Diag,S_W_B,S_D_B,N,m,(int)i1[0],(int)i2[0]);
}

