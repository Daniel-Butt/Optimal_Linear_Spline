#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm> 
#include <string>
#include "NumericMethods.h"


class OLS {

private:

    inline static double tol;

    static std::array<double, 2> fitLine(const double x1, const double x2, const double y1, const double y2) {

        //y = mx + b
        const double m = (y2 - y1) / (x2 - x1);
        const double b = y1 - x1 * m;

        return { m, b };
    } 

    static std::array<double, 2> splineWorstCaseError(std::function<double(double)> func, std::vector<double> I) {

        if (fabs(I[1] - I[0]) < 1e-7) return { 0.0, 0.0 };

        const auto& [m, b] = fitLine(I[0], I[1], func(I[0]), func(I[1]));

        //negate function to make maximum the minimum
        auto min_func = [&](double x) { return -abs(func(x) - (m * x + b)); };

        //negate output to make minimum the original maximum
        return NumericMethods::minimize(min_func, I, tol);
    }

    inline static std::function<double(double)> kmtempFunc;
    inline static double kmtempLbound;
    inline static double kmtempEps;

    static double knotRootFunc(double x) {
        return -splineWorstCaseError(kmtempFunc, { kmtempLbound, x })[1] - kmtempEps;
    }

    static double findNextKnot(std::function<double(double)> func, std::vector<double> I) {

        kmtempLbound = I[0];

        return NumericMethods::findRoot(knotRootFunc, I, tol);
    }

    static std::vector<double> olsMinmax(std::function<double(double)> func, std::vector<double> I, double eps) {

        kmtempFunc = func;
        kmtempEps = eps;

        double currentKnot = I[0];

        std::vector<double> knots = { currentKnot };

        while (I[1] - currentKnot > tol) {

            currentKnot = findNextKnot(func, { currentKnot, I[1] });
            knots.push_back(currentKnot);

        }

        return knots;
    }

    static std::pair<double, std::vector<double>> olsRefinedMinmax(std::function<double(double)> func, std::vector<double> I, double eps) {

        //binary search to refine epsilon (max error) value
        double lEps = 0.0;
        double uEps = eps;
        double currentEps = eps;
        double lastEps = 0.0;
        int minKnots = (int)olsMinmax(func, I, eps).size();

        while (abs(currentEps - lastEps) > tol) {

            lastEps = currentEps;
            currentEps = (lEps + uEps) / 2.0;
            const int num_knots = (int)olsMinmax(func, I, currentEps).size();

            if (num_knots > minKnots) {
                lEps = currentEps;
            }
            else {
                uEps = currentEps;
            }
        }

        return { uEps, olsMinmax(func, I, uEps) };
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

        LineNode operator[](const int i) {
            return table[i];
        }

        double operator()(const double x) {

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

        int size() {
            return table.size();
        }

        constexpr std::vector<LineNode>::iterator begin() {
            return table.begin();
        }

        constexpr std::vector<LineNode>::iterator end() {
            return table.end();
        }

        std::string toString() {
            std::string str = "";

            str += "Optimal Table{\n";
            for (const auto& line : table) {
                str += "   Line: {\n      min: " + std::to_string(line.minX) + "\n      max: " + std::to_string(line.maxX) + "\n      m: " + std::to_string(line.m) + "\n      b: " + std::to_string(line.b) +"\n      " + std::to_string(line.m) + "x + " + std::to_string(line.b) + "\n   }\n";
            }
            str += "}\n";

            return str;
        }
        
    };

    static std::pair<double, SplineTable> fitSpline(std::function<double(double)> func, std::vector<double> I, double eps = 0.01, const double t = 1e-7, const olsMethod method = refined, const bool adjustableLines = false) {

        tol = t;
        
        eps = adjustableLines ? 2.0 * eps : eps;

        std::pair<double, std::vector<double>> result;

        switch (method) {

            case standard:
                result = { eps, olsMinmax(func, I, eps) };
                break;

            case refined:
            default:
                result = olsRefinedMinmax(func, I, eps);
                break;
        }


        return { result.first, constructSplineTable(func, result.second, adjustableLines) };

    }

    static SplineTable constructSplineTable(std::function<double(double)> func, std::vector<double>& knots, const bool adjustableLines) {

        const int N = (int)knots.size() - 1;

        std::vector<LineNode> lines;
        lines.reserve(N);

        for (int i = 0; i < N; i++) {

            const auto& [m, b] = adjustableLines ? fitOptimalLine(func, { knots[i], knots[i + 1] }) : fitLine(knots[i], knots[i + 1], func(knots[i]), func(knots[i + 1]));

            lines.emplace_back(knots[i], knots[i + 1], m, b);
        }

        return SplineTable(lines);
    }

    static std::array<double, 2> fitOptimalLine(std::function<double(double)> func, std::vector<double> I) {

        const double x1 = I[0];
        const double x2 = I[1];
        const double y1 = func(x1);
        const double y2 = func(x2);

        const double m = (y2 - y1) / (x2 - x1);
        double b = y1 - x1 * m;

        const double xc = splineWorstCaseError(func, I)[0];

        const double h = 0.5 * (func(xc) - m * xc - b);

        return { m, b + h };
    }

    /*
        Mat3 A = { { {1.0, I[0], -1.0 },
                    {1.0,  0.0,  1.0 },
                    {1.0, I[1], -1.0 }  } };

        Vec3 b = { func(I[0]), 0.0, func(I[1]) };
        Vec3 c;
        double extreme = (I[0] + I[1]) / 2.0;
        double meanErr = std::numeric_limits<double>::max();

        auto R = [&](double x) { return func(x) - (c[1] * x + c[0]); };
        auto nabsR = [&](double x) { return -abs(R(x)); };

        for (int i = 0; i < 100; i++) {

            A[1][1] = extreme;
            b[1] = func(extreme);
            c = NumericMethods::linSolve3(A, b);

            std::vector<double> interval = { NumericMethods::findRoot(R, {I[0], extreme}), NumericMethods::findRoot(R, {extreme, I[1]}) };

            extreme = NumericMethods::minimize(nabsR, interval)[0];

            const Vec3 errs = { abs(R(I[0])), abs(R(extreme)), abs(R(I[1])) };

            meanErr = (errs[0] + errs[1] + errs[2]) / 3.0;

            const Vec3 meanDiffs = { abs(errs[0] - meanErr), abs(errs[1] - meanErr), abs(errs[2] - meanErr) };

            if (meanDiffs[0] < 1e-7 && meanDiffs[1] < 1e-7 && meanDiffs[2] < 1e-7) {
                break;
            }
        }

        return { c[1], c[0] };
    */
};
