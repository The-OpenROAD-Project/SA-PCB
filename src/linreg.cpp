///////////////////////////////////////////////////////////////////////////////
// Authors: Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng, Steven Swanson).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include "linreg.h"

inline static REAL sqr(REAL x) {
    return x*x;
}


int linreg(int n, const std::vector<double> x, const std::vector<double> y, double& m, double& b, double& r){
    REAL   sumx = 0.0;                      /* sum of x     */
    REAL   sumx2 = 0.0;                     /* sum of x**2  */
    REAL   sumxy = 0.0;                     /* sum of x * y */
    REAL   sumy = 0.0;                      /* sum of y     */
    REAL   sumy2 = 0.0;                     /* sum of y**2  */

    for (int i=0;i<n;i++){ 
        sumx  += x[i];       
        sumx2 += sqr(x[i]);  
        sumxy += x[i] * y[i];
        sumy  += y[i];      
        sumy2 += sqr(y[i]); 
    } 

    REAL denom = (n * sumx2 - sqr(sumx));
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        m = 0;
        b = 0;
        if (r) r = 0;
            return 1;
    }

    m = (n * sumxy  -  sumx * sumy) / denom;
    b = (sumy * sumx2  -  sumx * sumxy) / denom;
    

    REAL rmse;
    if(r!=NULL) {
        for (int i=0;i<n;i++){ 
            REAL lab = m * x[i] + b;
            rmse += sqrt(y[i] - lab);
        } 
        rmse = rmse / n;
        rmse = sqrt(rmse);
        r = rmse;
    }

    return 0; 
}