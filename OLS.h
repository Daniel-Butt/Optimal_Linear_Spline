#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm> 
#include "NumericOptimize.h"


class OLS {

private:

    static std::array<double, 2> fitLine(const double x1, const double x2, const double y1, const double y2) {

        //y = mx + b
        const double m = (y2 - y1) / (x2 - x1);
        const double b = y1 - x1 * m;

        return { m, b };
    }

    static double splineWorstCaseError(std::function<double(double)> func, std::vector<double> I) {

        if (fabs(I[1] - I[0]) < 1e-7) return 0;

        const auto& [m, b] = fitLine(I[0], I[1], func(I[0]), func(I[1]));

        //negate function to make maximum the minimum
        auto min_func = [&](double x) { return -abs(func(x) - (m * x + b)); };

        //negate output to make minimum the original maximum
        return -NumericOptimize::minimize(min_func, I)[1];
    }

    inline static std::function<double(double)> kmtempFunc;
    inline static double kmtempLbound;
    inline static double kmtempEps;

    static double knotRootFunc(double x) {
        return splineWorstCaseError(kmtempFunc, { kmtempLbound, x }) - kmtempEps;
    }

    static double findNextKnot(std::function<double(double)> func, std::vector<double> I, double eps, double tol = 1e-7) {

        kmtempLbound = I[0];

        return NumericOptimize::findRoot(knotRootFunc, I, tol);
    }

    static std::vector<double> olsMinmax(std::function<double(double)> func, std::vector<double> I, double eps, double tol) {

        kmtempFunc = func;
        kmtempEps = eps;

        double currentKnot = I[0];

        std::vector<double> knots = { currentKnot };

        while (I[1] - currentKnot > tol) {

            currentKnot = findNextKnot(func, { currentKnot, I[1] }, eps, tol);
            knots.push_back(currentKnot);

        }

        return knots;
    }

    static std::pair<double, std::vector<double>> olsRefinedMinmax(std::function<double(double)> func, std::vector<double> I, double eps, double tol) {

        //binary search to refine epsilon (max error) value
        double lEps = 0.0;
        double uEps = eps;
        double currentEps = eps;
        double lastEps = 0.0;
        int minKnots = (int)olsMinmax(func, I, eps, tol).size();

        while (abs(currentEps - lastEps) > tol) {

            lastEps = currentEps;
            currentEps = (lEps + uEps) / 2.0;
            const int num_knots = (int)olsMinmax(func, I, currentEps, tol).size();

            if (num_knots > minKnots) {
                lEps = currentEps;
            }
            else {
                uEps = currentEps;
            }
        }

        return { uEps, olsMinmax(func, I, uEps, tol) };
    }


public:

    enum olsMethod {
        standard,
        refined
    };

    struct LineNode
    {
        double minX;
        double maxX;
        double m;
        double b;
    };

    struct SplineTable {
        std::vector<LineNode> table;

        double operator[](const double x) {

            //binary search for correct line
            int lBound = 0;
            int uBound = (int)table.size();

            while (true) {
                const int currentIdx = (uBound + lBound) >> 1;
                const LineNode& node = table[currentIdx];

                if (x < node.minX) {
                    uBound = currentIdx;
                }
                else if (x > node.maxX) {
                    lBound = currentIdx;
                }
                else {
                    //compute linear interpolation
                    return node.m * x + node.b;
                }
            }
        }
    };

    static std::pair<double, std::vector<double>> fitSpline(std::function<double(double)> func, std::vector<double> I, const olsMethod method = standard, const double eps = 0.01, const double tol = 1e-7) {

        switch (method) {

            case standard:
                return { eps, olsMinmax(func, I, eps, tol) };

            case refined:
                return olsRefinedMinmax(func, I, eps, tol);

            default:
                return { eps, olsMinmax(func, I, eps, tol) };
        }
    }

    static SplineTable constructSplineTable(std::function<double(double)> func, std::vector<double>& knots) {

        const int N = (int)knots.size() - 1;

        std::vector<LineNode> lines;
        lines.reserve(N);

        for (int i = 0; i < N; i++) {

            const auto& [m, b] = fitLine(knots[i], knots[i + 1], func(knots[i]), func(knots[i + 1]));

            lines.emplace_back(knots[i], knots[i + 1], m, b);
        }

        return SplineTable(lines);
    }

};
