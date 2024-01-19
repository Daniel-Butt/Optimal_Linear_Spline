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

double log2_(double x) {

    return std::log2(x);

}


double recp(double x) {
    return 1.0 / x;
}

double foo(double x) {
    return 4.3 * pow(x, 3) - 2.3 * pow(x, 2) - 6.6 * x + 2.6; // [0.379, 1.343] 
}

double SIN(double x) {
    return sin(x);
}

double burgersRottTan(double x) {
    const double L = 8.78359490542964;
    const double a = 1.2564312086261704;
    const double pi = 3.141592653589793;


    return x < 1e-7 ? 0.0 : (L / (2.0 * pi * x)) * (1.0 - exp(-a * x * x));
}
 
int main(){

    auto result = OLS::fitSpline(recp, { 1.0, 2.0 }, 0.001, 1e-10, OLS::refined, true);

    std::cout << result.second.toString() << std::endl;

    //std::cout << result.first << ", " << result.second.size() << std::endl;

    /*for (const auto& x : result.second)
        std::cout << x << std::endl;*/

    //auto funcApprox = OLS::constructSplineTable(log2_, result.second);

    //const double x = 1.467;
 
    //printf("Real: %.16f,  Approx: %.16f\n\n", log2_(x), funcApprox(x));

    ////std::cout << sinApprox.toString();


    /*auto result = OLS::fitOptimalLine(recp, { 1.0, 2.0 });

    for (const auto& x : result)
        std::cout << x << std::endl;*/
}

