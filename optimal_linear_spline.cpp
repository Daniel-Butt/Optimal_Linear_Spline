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

double SIN(double x) {
    return sin(x);
}
 
int main(){

    auto result = OLS::fitSpline(SIN, { 0.0, (3.141592653589793) / 2.0 }, OLS::refined, 0.001, 1e-14).second;

    /*std::cout << result.first << ", " << result.second.size() << std::endl;

    for (const auto& x : result.second)
        std::cout << x << std::endl;*/

    auto sinApprox = OLS::constructSplineTable(SIN, result);

    const double x = 0.562;
 
    printf("Real: %.5f,  Approx: %.5f\n\n", sin(x), sinApprox(x));

    std::cout << sinApprox.toString();
}

