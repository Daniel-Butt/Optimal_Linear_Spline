#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include "NumericMethods.h"

class ProfileSpline
{

private:
	double eps;

	double binarySearch(std::function<double(double)> func, std::array<double, 2> I, double tol = 1e-7) {
		
		double& left = I[0];
		double& right = I[1];
		double mid = (right + left) * 0.5;

		while (right - left > tol) {

			if (func(mid) > 0.0) {
				left = mid;
			}
			else {
				right = mid;
			}

			mid = (right + left) * 0.5;
		}
		return mid;
	}

	std::array<double, 2> line(double x1, double y1, double x2, double y2) {
		double m = (y2 - y1) / (x2 - x1);
		double b = y1 - m * x1;

		return { m, b };
	}

	std::array<double, 2> knotSearch(std::array<double, 2> startPoint, double end, double tol = 1e-7) {
		
		double dir = 1.0;

		auto errorFunc = [&](double x) {

			auto [ma, ba] = line(startPoint[0], startPoint[1], x, func(x) + eps);
			auto aboveErrorFunc = [&](double x1) {
				return -fabs(func(x1) - (ma * x1 + ba));
			};

			auto [mb, bb] = line(startPoint[0], startPoint[1], x, func(x) - eps);
			auto belowErrorFunc = [&](double x1) {
				return -fabs(func(x1) - (mb * x1 + bb));
			};

			const double aboveError = -NumericMethods::globalMinimize(aboveErrorFunc, { startPoint[0], x }, 0.01)[1];
			const double belowError = -NumericMethods::globalMinimize(belowErrorFunc, { startPoint[0], x }, 0.01)[1];

			dir = aboveError < belowError ? 1.0 : -1.0;

			return std::min(aboveError, belowError) - eps;
		};

		const double x = NumericMethods::findRoot(errorFunc, { startPoint[0], end }, tol);
		
		return { x, func(x) + eps * dir };
	}

	double worstErrorAtEnd(std::array<double, 2> startPoint, double end, double tol = 1e-7) {
		auto [m, b] = line(startPoint[0], startPoint[1], end, func(end));

		auto errorFunc = [&](double x) {
			return -fabs(func(x) - (m * x + b));
		};

		const double error = -NumericMethods::globalMinimize(errorFunc, { startPoint[0], end }, 0.01)[1];

		return error;
	}

	std::vector<std::array<double, 2>> fitSpline(std::array<double, 2> startPoint, double end) {
		std::vector<std::array<double, 2>> points;

		std::array<double, 2> p = startPoint;
		points.push_back(p);

		while (worstErrorAtEnd(p, end) > eps) {
			p = knotSearch(p, end);
			points.push_back(p);
		}

		return points;
	}

	std::vector<std::array<double, 2>> fitRefinedSpline(std::array<double, 2> startPoint, double end) {
		auto pointsFunc = [&](double epsilon) {
			eps = epsilon;
			std::vector<std::array<double, 2>> points = fitSpline(startPoint, end);
			return points.size();
		};

		int numPointsRequired = pointsFunc(e);

		auto pointsErrorFunc = [&](double epsilon) {
			double p = pointsFunc(epsilon);
			return p - numPointsRequired;
		};

		eps = binarySearch(pointsErrorFunc, { 0.0, e }, 1e-7) + 1e-7;

		std::vector<std::array<double, 2>> points = fitSpline(startPoint, end);

		points.push_back({ end, func(end) });

		return points;
	}

public:

	std::vector<std::array<double, 2>> points;
	std::function<double(double)> func;
	double e;
	
	ProfileSpline(std::function<double(double)> f, double e_) : func(f), e(e_) {

		auto innerPoints = fitRefinedSpline({ 0.0, 0.0 }, 1.0);
		auto outerPoints = fitRefinedSpline({ 1.0, 1.0 }, 10.0);

		points.insert(points.end(), innerPoints.begin(), innerPoints.end());
		points.insert(points.end(), outerPoints.begin()+1, outerPoints.end());
		
		for (const auto& p : points) {
			std::cout << p[0] << ", " << p[1] << std::endl;
		}
	}

};
