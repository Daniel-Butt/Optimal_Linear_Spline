#include <iostream>
#include "ProfileSpline.h"
#include "SullivanLookUp.h"
#include <chrono>

double f(double x){
	//return x < 1 ? pow(x, 0.5) : pow(x, -0.5);
	return 2 * x / (x * x + 1);
	//constexpr double G = 8.783595;
	//constexpr double a = 1.2564312;
	//constexpr double PI = 3.14159265358979323846;

	//return G * (1.0 - exp(-a * x * x)) / (2.0 * PI * x);
}

double g(double x) {
	return sin(x) + sin(10.0 / 3.0 * x);
}

double phi;

double rankine(double x) {

	return x < 1 ? x : 1.0/x;
}


int main(){

	ProfileSpline ps(f, 0.01);

	for (int i = 1; i < ps.points.size(); i++) {
		const double x1 = ps.points[i - 1][0];
		const double x2 = ps.points[i][0];
		const double y1 = ps.points[i - 1][1];
		const double y2 = ps.points[i][1];

		const double m = (y2 - y1) / (x2 - x1);
		const double b = y1 - m * x1;

		printf("{ %f, %f, %f, %f, %f },\n", x2, m, b, m, b);
	}
}

