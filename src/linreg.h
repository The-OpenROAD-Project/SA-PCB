#include <stdlib.h>
#include <math.h>                           /* math functions */

//#define REAL float
#define REAL double

int linreg(int n, const std::vector<double> x[], const std::vector<double> y[], double& m, double& b, double& r);


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