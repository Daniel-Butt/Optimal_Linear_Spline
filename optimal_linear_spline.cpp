#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <functional>
#include "OLS.h"


//test functions
double f(double x) {
    return std::pow(x, 0.6);
}

double foo(double x) {
    return 4.3 * pow(x, 3) - 2.3 * pow(x, 2) - 6.6 * x + 2.6; // [0.379, 1.343] 
}

 
int main(){

    auto result = OLS::fitSpline(f, { 0.0, 1.0 }, OLS::refined, 0.001, 1e-7).second;

    auto table = OLS::constructSplineTable(f, result);

    const double x = 0.562;
 
    printf("Real: %.4f,  Approx: %.4f\n", f(x), OLS::compute(table, x));


}

