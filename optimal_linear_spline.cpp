#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <functional>
#include "NumericOptimize.h"


//test functions
double f(double x) {
    return std::pow(x, 0.6);
}

double foo(double x) {
    return 4.3 * pow(x, 3) - 2.3 * pow(x, 2) - 6.6 * x + 2.6; // [0.379, 1.343] 
}

std::array<double, 2> fitLine(const double x1, const double x2, const double y1, const double y2) {

    //y = mx + b
    const double m = (y2 - y1) / (x2 - x1);
    const double b = y1 - x1 * m;

    return { m, b };
}

double spline_worst_case_error(std::function<double(double)> func, std::vector<double> I) {

    if (fabs(I[1] - I[0]) < 1e-7) return 0;

    const auto& [m, b] = fitLine(I[0], I[1], func(I[0]), func(I[1]));

    //negate function to make maximum the minimum
    auto min_func = [&](double x) { return -abs( func(x) - (m * x + b) ); };

    //negate output to make minimum the original maximum
    return -NumericOptimize::minimize(min_func, I)[1];
}


static std::function<double(double)> kmtemp_func;
static double kmtemp_lbound;
static double kmtemp_eps;

double knot_root_func(double x) {
    return spline_worst_case_error(kmtemp_func, { kmtemp_lbound, x }) - kmtemp_eps;
}

double find_next_knot(std::function<double(double)> func, std::vector<double> I, double eps, double tol = 1e-7) {

    kmtemp_lbound = I[0];

    return NumericOptimize::findRoot(knot_root_func, I, tol);

    /*double l_bound = I[0];
    double u_bound = I[1];
    double error = spline_worst_case_error(func, I);
    double knot = u_bound;

    if (error < eps) {
        return knot;
    }

    while (std::abs(error - eps) > tol) {

        knot = (u_bound + l_bound) / 2.0;
        error = spline_worst_case_error(func, { I[0], knot });

        if (error > eps) {
            u_bound = knot;
        }
        else {
            l_bound = knot;
        }
    }

    return knot;*/
}

std::vector<double> ols_minmax(std::function<double(double)> func, std::vector<double> I, double eps, double tol) {

    kmtemp_func = func;
    kmtemp_eps = eps;

    double current_knot = I[0];

    std::vector<double> knots = { current_knot };

    while (I[1] - current_knot > tol) {

        current_knot = find_next_knot(func, { current_knot, I[1] }, eps, tol);
        knots.push_back(current_knot);

    }

    return knots;
}

std::pair<double, std::vector<double>> ols_refined_minmax(std::function<double(double)> func, std::vector<double> I, double eps, double tol) {

    //binary search to refine epsilon (max error) value
    double l_eps = 0.0;
    double u_eps = eps;
    double current_eps = eps;
    double last_eps = 0.0;
    int min_knots = ols_minmax(func, I, eps, tol).size();

    while (abs(current_eps - last_eps) > tol) {

        last_eps = current_eps;
        current_eps = (l_eps + u_eps) / 2.0;
        const int num_knots = ols_minmax(func, I, current_eps, tol).size();

        if (num_knots > min_knots) {
            l_eps = current_eps;
        }
        else {
            u_eps = current_eps;
        }
    }

    return { u_eps, ols_minmax(func, I, u_eps, tol) };
}


enum olsMethod {
    standard,
    refined
};

std::pair<double, std::vector<double>> optimal_linear_spline(std::function<double(double)> func, std::vector<double> I, const olsMethod method = standard, const double eps = 0.01, const double tol = 1e-7) {

    switch (method) {

        case standard:
            return { eps, ols_minmax(func, I, eps, tol) };

        case refined:
            return ols_refined_minmax(func, I, eps, tol);
    }
}


int main(){

    const auto result = optimal_linear_spline(f, { 0.0, 1.0 }, refined, 0.05);

    printf("%.7f\n\n", result.first);

    printf("[\n");

    for (const auto& x : result.second) {
        printf("    %.7f,\n", x);
    }

    printf("]\n");



}

