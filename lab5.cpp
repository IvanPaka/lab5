#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void spline(double *fx, int n, double *c, double *b, double *d){
    int m = n - 2;
    double a[m][m+1];
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            if(i == j && j-1 >= 0) {
                a[i][j] = 4.;
                a[i][j-1] = 1.;
                a[i][j+1] = 1.;
            }else if(i == j) {
                a[i][j] = 4.;
                a[i][j+1] = 1.;
            }else if(j > i+1 || j < i-1)
                a[i][j] = 0.;
        }
        a[i][m] = 3.*(fx[i+2] - 2.*fx[i+1] + fx[i]);
    }

    for(int k = 1; k < m; k++){
        for(int j = k; j < m; j++){
            double h = a[j][k-1]/a[k-1][k-1];
            for(int i = 0; i < m+1; i++)
                a[j][i] -= h*a[k-1][i];
        }
        for(int i = m-1; i >= 0; i--){
            c[i+1] = a[i][m]/a[i][i];
            for(int l = m-1; l > i; l--)
                c[i+1] -= a[i][l]*c[l+1]/a[i][i];
        }
    }

    for(int i = 0; i < m+1; i++) {
        b[i] = fx[i+1] - fx[i] - (c[i+1] + 2*c[i])/3.;
        d[i] = (c[i+1] - c[i])/3.;
        if(i == m){
            b[i] = fx[i+1] - fx[i] - 2*c[i]/3.;
            d[i] = -c[i]/3.;
        }
    }
}



double trapez(double *f, int size, double h) {
    double result = h*0.5*(f[0] + f[size-1]);
    int i = 0;
    for(int i = 1; i < size - 1; i++){
        result += h*f[i];
        i++;
    }
    return result;
}


int main(){
    int n = 9;
    double h = 0.25;
    double fx[] = {-0.333, 0, -0.125, -0.056, 0, 0.046, 0.083, 0.115, 0.143};
    double c[n-1], b[n-1], d[n-1];
    spline(fx, n, c, b, d);


    double p[(n-1)*20 + 1];
    int size = (n-1)*20 + 1;
    for(int i = 0; i < n-1; i++){
            for(int j = 0; j < 20; j++){
                double xx = j*0.1;
                p[i*20+j] = fx[i] + b[i]*xx + c[i]*xx*xx + d[i]*xx*xx*xx;
            }
        }
    p[(n-1)*20] = fx[n-1];


    double trapez_result = trapez(p, size, h/20.);

    double p2[(n-1)*10+1];
    size = (n-1)*10 + 1;
    for(int i = 0; i < n-1; i++){
            for(int j = 0; j < 10; j++){
                double xx = j*0.2;
                p2[i*10+j] = fx[i] + b[i]*xx + c[i]*xx*xx + d[i]*xx*xx*xx;
            }
        }
    p2[(n-1)*10 + 1] = fx[n-1];

    double trapez_result2 = trapez(p2, size, h/10.);

    ofstream outf;
    outf.open("ans1.dat", ios_base::out);
    outf << trapez_result << " " << abs(trapez_result - trapez_result2)/15.;
    outf.close();

    return 0;
}
