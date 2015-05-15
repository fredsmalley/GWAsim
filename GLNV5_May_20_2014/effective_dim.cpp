//
//  effective_dim.cpp
//  gln
//
//  Created by Joe Song on 3/10/13.
//
//
#include <vector>
using namespace std;

#include "nr3.h"
#include "eigen_sym.h"

size_t cal_eff_dim(const vector< vector<int> > & X)
// compute the effective dimension of X
// X: NxD matrix. N: sample size. D: dimension
{
    size_t D_eff;
    
    size_t N = X.size();
    size_t D;
    
    if(N == 0) {
        
        D_eff = 0;
        return D_eff;
        
    } else if (X[0].size() > 1000) {
        
        D = X[0].size();
        D_eff = N-1 < D ? N-1 : D;
        return D_eff;

    } else {
        
        D = X[0].size();
        D_eff = D;
    }
    
    vector<double> mean(D, 0), sd(D, 0);
    for(size_t i=0; i<D; i++) {
        for(size_t n=0; n<N; n++) {
            mean[i] += X[n][i];
            sd[i] += X[n][i] * X[n][i];
        }
        mean[i] /= N;
        sd[i] = sqrt( sd[i] / N - mean[i] * mean[i] );
    }
    
    MatDoub C(D, D);
    
    // Compute correlation matrix C for X
    for(size_t i=0; i<D; i++) {
        for(size_t j=0; j<=i; j++) {
            C[i][j] = 0.0;
            for(size_t n=0; n<N; n++) {
                // C[i][j] += (X[n][i] * X[n][j] - mean[i] * mean[j]) / (sd[i] * sd[j]);
                C[i][j] += X[n][i] * X[n][j];
            }
            C[i][j] = ( C[i][j] / N - mean[i] * mean[j] ) / (sd[i] * sd[j]);
        }
    }
    
    for(size_t i=0; i<D; i++) {
        for(size_t j=i+1; j < D; j++) {
            C[i][j] = C[j][i];
        }
    }
    
    // Compute eigen-value of C
    Symmeig s(C, false);
    
    // Compute the effective dimesion based on C
    double sum = 0.0;
    for(size_t i=0; i<D; i++) {
        sum += abs(s.d[i]);
        if (sum > D-1) {
            D_eff = i+1;
            break;
        }
    }
    
    return D_eff;
}