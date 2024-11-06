#pragma once

#include <cmath>
#include <functional>
#include <array>
#include <vector>
#include <complex>

typedef std::complex<double> Complex;
typedef std::vector<std::vector<double>> Mat;
typedef std::vector<double> Vec;
typedef std::vector<Complex> CVec;
typedef std::vector<int> IntVec;
typedef std::array<std::array<double, 3>, 3> Mat3;
typedef std::array<double, 3> Vec3;
typedef std::array<Complex, 3> CVec3;
typedef std::array<std::array<double, 2>, 2> Mat2;
typedef std::array<double, 2> Vec2;
typedef std::array<Complex, 2> CVec2;


class NumericMethods {

private:
    static double goldMinimize(std::function<double(double)> f, std::array<double, 2>& I, const double tol = 1e-7) {

        const double phi = (sqrt(5.0) + 1.0) / 2.0; // Golden ratio
        double a = I[0];
        double b = I[1];

        while (abs(b - a) > tol) {

            double c = b - (b - a) / phi;
            double d = a + (b - a) / phi;

            if (f(c) < f(d)) {
                b = d;
            }
            else {
                a = c;
            }
        }
        return (b + a) / 2.0;
    }

    static double brentMinimize(std::function<double(double)> f, std::array<double, 2>& I, const double t = 1e-7) {
        double a = I[0];
        double b = I[1];

        double c, d, e, eps, fu, fv, fw, fx, m, p, q, r, sa, sb, t2, tol, u, v, w, x;

        // C is the square of the inverse of the golden ratio.
        c = 0.5 * (3.0 - sqrt(5.0));

        eps = sqrt(2.220446049250313E-016);

        sa = a;
        sb = b;
        x = sa + c * (b - a);
        w = x;
        v = w;
        d = 0.0;
        e = 0.0;
        fx = f(x);
        fw = fx;
        fv = fw;

        while (true) {
            m = 0.5 * (sa + sb);
            tol = eps * fabs(x) + t;
            t2 = 2.0 * tol;

            // Check the stopping criterion.
            if (fabs(x - m) <= t2 - 0.5 * (sb - sa))
            {
                break;
            }

            // Fit a parabola.
            r = 0.0;
            q = r;
            p = q;

            if (tol < fabs(e)) {
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);

                if (0.0 < q) {
                    p = -p;
                }

                q = fabs(q);
                r = e;
                e = d;
            }

            if (fabs(p) < fabs(0.5 * q * r) && q * (sa - x) < p && p < q * (sb - x)) {

                // Take the parabolic interpolation step.
                d = p / q;
                u = x + d;

                // F must not be evaluated too close to A or B.
                if ((u - sa) < t2 || (sb - u) < t2) {
                    if (x < m) {
                        d = tol;
                    }
                    else {
                        d = -tol;
                    }
                }
            }

            // A golden-section step.
            else {
                if (x < m) {
                    e = sb - x;
                }
                else {
                    e = sa - x;
                }
                d = c * e;
            }

            // F must not be evaluated too close to X.
            if (tol <= fabs(d)) {
                u = x + d;
            }
            else if (0.0 < d) {
                u = x + tol;
            }
            else {
                u = x - tol;
            }

            fu = f(u);

            // Update A, B, V, W, and X.
            if (fu <= fx) {
                if (u < x) {
                    sb = x;
                }
                else {
                    sa = x;
                }

                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            }
            else {
                if (u < x) {
                    sa = u;
                }
                else {
                    sb = u;
                }

                if (fu <= fw || w == x) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }

        return x;
    }

    static double brentDekkerRoot(std::function<double(double)> func, std::array<double, 2>& I, const double tol = 1e-7) {
        int iter;

        double a = I[0], b = I[1], c = I[1], d = 0.0, e = 0.0, min1, min2;
        double fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;

        //not bracketed
        if (fa * fb > 0.0) return b;

        fc = fb;

        for (iter = 1; iter <= 1000; iter++) {
            // if sign(fb) = sign(fc)
            if (fb * fc > 0.0) {
                c = a;
                fc = fa;
                e = d = b - a;
            }
            if (fabs(fc) < fabs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            tol1 = 2.0 * sqrt(2.220446049250313E-016) * fabs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);

            if (fabs(xm) <= tol1 || fb == 0.0) return b;

            if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                s = fb / fa;

                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) q = -q;

                p = fabs(p);

                min1 = 3.0 * xm * q - fabs(tol1 * q);
                min2 = fabs(e * q);

                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d;
                    d = p / q;
                }
                else {
                    d = xm;
                    e = d;
                }
            }
            else {
                d = xm;
                e = d;

            }
            a = b;
            fa = fb;

            if (fabs(d) > tol1)
                b += d;
            else
                b += copysign(tol1, xm);

            fb = func(b);
        }
        printf("Maximum number of iterations exceeded in brent root finding\n");
        return 0.0;
    }

    static inline double det3(const Mat3& mat) {

        return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
               mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
               mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    }

    static inline double det2(const Mat2& mat) {
        return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    }

    struct Adapt {
        const double toler;
        bool terminate;
        bool outOfTolerance;
    };

    static double adaptlob(std::function<double(double)> func, const double a, const double b, const double fa, const double fb, const double is, Adapt& adapt)
    {
        constexpr double alpha = 0.816496580927726;
        constexpr double beta = 0.4472135954999579;

        const double m = 0.5 * (a + b);
        const double h = 0.5 * (b - a);

        const double mll = m - alpha * h;
        const double ml = m - beta * h;
        const double mr = m + beta * h;
        const double mrr = m + alpha * h;

        const double fmll = func(mll);
        const double fml = func(ml);
        const double fm = func(m);
        const double fmr = func(mr);
        const double fmrr = func(mrr);

        const double i2 = h / 6.0 * (fa + fb + 5.0 * (fml + fmr));
        const double i1 = h / 1470.0 * (77.0 * (fa + fb) + 432.0 * (fmll + fmrr) + 625.0 * (fml + fmr) + 672.0 * fm);

        if (abs(i1 - i2) <= adapt.toler * is || mll <= a || b <= mrr) {

            if ((mll <= a || b <= mrr) && terminate) {
                adapt.outOfTolerance = true;
                adapt.terminate = false;
            }

            return i1;
        }

        else {
            return adaptlob(func, a, mll, fa, fmll, is, adapt) +
                   adaptlob(func, mll, ml, fmll, fml, is, adapt) +
                   adaptlob(func, ml, m, fml, fm, is, adapt) +
                   adaptlob(func, m, mr, fm, fmr, is, adapt) +
                   adaptlob(func, mr, mrr, fmr, fmrr, is, adapt) +
                   adaptlob(func, mrr, b, fmrr, fb, is, adapt);
        }
    }

    static Vec chebNodes(Vec2& I, const int n) {
        constexpr double PI = 3.141592653589793;

        const double k1 = 0.5 * (I[1] + I[0]);
        const double k2 = 0.5 * (I[1] - I[0]);
        const double k3 = PI / (2 * n);

        Vec nodes;
        nodes.reserve(n + 2);

        nodes.push_back(I[0]);

        for (int i = n; i >= 1; i--) {
            nodes.push_back(k1 + k2 * cos(k3 * (2 * i - 1)));
        }

        nodes.push_back(I[1]);

        return nodes;
    }

    static Mat constructMinimaxMat(Vec2& I, const int n) {

        Mat A;
        A.reserve(n + 2);

        Vec nodes = chebNodes(I, n);

        for (int i = 0; i < n + 2; i++) {

            A.emplace_back(Vec());
            A[i].resize(n + 2);
            A[i][0] = 1.0;
            A[i][1] = nodes[i];
            A[i][n + 1] = -pow(-1, i);
        }

        return A;
    }


public:

    static double sgn(double x) {
        return x > 0 ? 1.0 : -1.0;
    }

    static double findRoot(std::function<double(double)> func, std::array<double, 2> I, const double tol = 1e-7) {
        return brentDekkerRoot(func, I, tol);
    }

    static Vec2 minimize(std::function<double(double)> func, std::array<double, 2> I, const double tol = 1e-7) {
        const double minX = brentMinimize(func, I, tol);

        return { minX, func(minX) };
    }

    static Vec2 maximize(std::function<double(double)> func, std::array<double, 2> I, const double tol = 1e-7) {

        const double maxX = brentMinimize([&](double x) { return -func(x); }, I, tol);

        return { maxX, func(maxX) };
    
    }

    /*static Vec2 bruteGlobalMinimum(std::function<double(double)> func, std::array<double, 2> I, const double tol = 1e-4) {
    
        //assumes function is continuous and multimodal

        double minVal = 0.0;
        double minX = 0.0;

        const double delta = (I[1] - I[0]) * tol;

        for (double x = I[0]; x <= I[1]; x += delta) {

            const double d = func(x);

            if (d < minVal) {
                minVal = d;
                minX = x;
            }
        }

        return NumericMethods::minimize([&](double x) { return func(x); }, { std::max(I[0], minX - delta), std::min(I[1], minX + delta)}, tol);
    }

    static Vec2 bruteGlobalMaximum(std::function<double(double)> func, std::array<double, 2> I, const double tol = 1e-4) {

        //assumes function is continuous and multimodal

        double maxVal = 0.0;
        double maxX = 0.0;

        const double delta = (I[1] - I[0]) * tol;

        for (double x = I[0]; x <= I[1]; x += delta) {

            const double d = func(x);

            if (d > maxVal) {
                maxVal = d;
                maxX = x;
            }
        }

        return NumericMethods::maximize([&](double x) { return func(x); }, { std::max(I[0], maxX - delta), std::min(I[1], maxX + delta) }, tol);
    }*/

    static std::vector<Vec2> localMinimums(std::function<double(double)> func, std::array<double, 2> I, const double mcr = 1e-2, const double tol = 1e-7) {
        
        std::vector<Vec2> localMinima;
        const int n = 4 * std::ceil((I[1] - I[0]) / mcr);

        std::vector<double> x, y;
        x.reserve(n);
        y.reserve(n);

        double step = (I[1] - I[0]) / (n - 1);

        for (int i = 0; i < n; ++i) {
            x.push_back(I[0] + i * step);
            y.push_back(func(x[i]));
        }

        std::vector<double> dy;
        dy.reserve(n - 1);

        for (int i = 0; i < n - 1; ++i) {
            dy.push_back(y[i + 1] - y[i]);
        }

        for (int i = 0; i < n - 2; ++i) {
            if (dy[i] < 0 && dy[i + 1] > 0) {
                std::array<double, 2> subI = { x[i], x[i + 2] };
                const double xmin = brentMinimize(func, subI, tol);
                const double fval = func(xmin);

                localMinima.push_back({ xmin, fval });
            }
        }
    
        return localMinima;
    }

    static Vec2 globalMinimize(std::function<double(double)> func, std::array<double, 2> I, const double mcr = 1e-2, const double tol = 1e-7) {
    
        double global_min = 1e99;
        double global_min_x = std::nan("");

        if (fabs(I[1] - I[0]) < mcr) {
			return { I[0], 0.0 };
		}

        const int n = std::max((int)(4 * (I[1] - I[0]) / mcr) + 1, 3);

        std::vector<double> x, y;
        x.reserve(n);
        y.reserve(n);

        double step = (I[1] - I[0]) / (n - 1);

        for (int i = 0; i < n; ++i) {
            x.push_back(I[0] + i * step);
            y.push_back(func(x[i]));
        }

        std::vector<double> dy;
		dy.reserve(n - 1);

        for (int i = 0; i < n - 1; ++i) {
            dy.push_back(y[i + 1] - y[i]);
        }

        for (int i = 0; i < n - 2; ++i) {
            if (dy[i] < 0 && dy[i + 1] > 0) {

                std::array<double, 2> subI = { x[i], x[i + 2] };
				const double xmin = brentMinimize(func, subI, tol);
                const double fval = func(xmin);

                if (fval < global_min) {
                    global_min = fval;
                    global_min_x = xmin;
                }
            }
        }
        
		return { global_min_x, global_min };
    }


    //Crammer's Rule
    static Vec3 linSolve3(Mat3& A, Vec3& b) {

        const double det_1 = 1.0 / det3(A);

        const Mat3& mi = {{ {b[0], A[0][1], A[0][2]},
                            {b[1], A[1][1], A[1][2]},
                            {b[2], A[2][1], A[2][2]}  }};

        const Mat3& mj = {{ {A[0][0], b[0], A[0][2]},
                            {A[1][0], b[1], A[1][2]},
                            {A[2][0], b[2], A[2][2]}  }};

        const Mat3& mk = {{ {A[0][0], A[0][1], b[0]},
                            {A[1][0], A[1][1], b[1]},
                            {A[2][0], A[2][1], b[2]}  }};


        return { det3(mi) * det_1, det3(mj) * det_1, det3(mk) * det_1 };

    }

    static Vec2 linSolve2(Mat2& A, Vec2& b) {

        const double det_1 = 1.0 / det2(A);

        Mat2 M0 = {{ { b[0], A[0][1] },
                     { b[1], A[1][1] }  }};

        Mat2 M1 = {{ { A[0][0], b[0] },
                     { A[1][0], b[1] }  }};

        return { det2(M0) * det_1, det2(M1) * det_1 };
    }

    //LU decomposition
    static Vec linSolve(Mat A, Vec b) {

        const int n = (int) b.size();

        const double TINY = 1.0e-40;
        int imax;
        double big, temp;

        Vec vv;
        vv.resize(n);

        IntVec indx;
        indx.resize(n);

        for (int i = 0; i < n; i++) {

            big = 0.0;

            for (int j = 0; j < n; j++) {

                if ((temp = abs(A[i][j])) > big) big = temp;
            }

            if (big == 0.0) return { nan("") }; //singular matrix
 
            vv[i] = 1.0 / big; 
        }

        for (int k = 0; k < n; k++) {

            big = 0.0;

            for (int i = k; i < n; i++) {
                temp = vv[i] * abs(A[i][k]);

                if (temp > big) {
                    big = temp;
                    imax = i;
                }
            }

            if (k != imax) {
                for (int j = 0; j < n; j++) {

                    temp = A[imax][j];
                    A[imax][j] = A[k][j];
                    A[k][j] = temp;

                }
                vv[imax] = vv[k];
            }

            indx[k] = imax;
            if (A[k][k] == 0.0) A[k][k] = TINY;

            for (int i = k + 1; i < n; i++) {
                temp = A[i][k] /= A[k][k];

                for (int j = k + 1; j < n; j++) {

                    A[i][j] -= temp * A[k][j];
                }
            }
        }

        int ii = 0, ip;
        double sum;
        Vec x;
        x.resize(n);

        for (int i = 0; i < n; i++) x[i] = b[i];

        for (int i = 0; i < n; i++) {

            ip = indx[i];
            sum = x[ip];
            x[ip] = x[i];

            if (ii != 0) {
                for (int j = ii - 1; j < i; j++) sum -= A[i][j] * x[j];
            }
            else if (sum != 0.0) {
                ii = i + 1;
            }

            x[i] = sum;
        }

        for (int i = n - 1; i >= 0; i--) {

            sum = x[i];

            for (int j = i + 1; j < n; j++) sum -= A[i][j] * x[j];

            x[i] = sum / A[i][i];
        }

        return x;
    }

    static double integrate(std::function<double(double)> func, const double a, const double b, const double tol = 1e-7)
    {
        constexpr double alpha = 0.816496580927726;
        constexpr double beta = 0.4472135954999579;
        constexpr double x1 = 0.942882415695480;
        constexpr double x2 = 0.641853342345781;
        constexpr double x3 = 0.236383199662150;
        constexpr double x[12] = { 0.0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1 };

        double y[13];

        const double m = 0.5 * (a + b);
        const double h = 0.5 * (b - a);

        const double fa = y[0] = func(a);
        const double fb = y[12] = func(b);

        for (int i = 1; i < 12; i++) {

            y[i] = func(m + x[i] * h);
        }

        const double i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));
        const double i1 = (h / 1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]);

        double is = h * (0.0158271919734802 * (y[0] + y[12]) + 0.0942738402188500 * (y[1] + y[11]) + 0.155071987336585 * (y[2] + y[10]) +
                    0.188821573960182 * (y[3] + y[9]) + 0.199773405226859 * (y[4] + y[8]) + 0.224926465333340 * (y[5] + y[7]) + 0.242611071901408 * y[6]);

        const double erri1 = abs(i1 - is);
        const double erri2 = abs(i2 - is);

        const double r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;
        const double toler = (r > 0.0 && r < 1.0) ? tol / r : tol;

        if (is == 0.0) is = b - a;

        is = abs(is);

        Adapt adapt = { toler, true, false };

        return adaptlob(func, a, b, fa, fb, is, adapt);
    }

    static double differentiate(std::function<double(double)> func, const double x, const double h)
    {
        constexpr int ntab = 10;
        constexpr double k1 = 0.7142857142857143, k2 = 1.96;
        constexpr double big = std::numeric_limits<double>::max();
        constexpr double safe = 2.0;

        double errt, fac;
        double ans = 0.0;

        Mat a;
        a.resize(ntab);

        for (int i = 0; i < ntab; i++) a[i].resize(ntab);


        double hh = h;

        a[0][0] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

        double err = big;

        for (int i = 1; i < ntab; i++) {

            hh *= k1;

            a[0][i] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

            fac = k2;

            for (int j = 1; j <= i; j++) {

                a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);

                fac = k2 * fac;

                errt = std::max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]));

                if (errt <= err) {
                    err = errt;
                    ans = a[j][i];
                }
            }
            if (abs(a[i][i] - a[i - 1][i - 1]) >= safe * err) break;
        }
        return ans;
    }

    static double differentiate(std::function<double(double)> func, const double x){
        return differentiate(func, x, 1e-3);
    }

    //coeffs in acending order x^0, x^1, ... x^n
    static double polyEval(Vec c, const double x) {
        const int n = c.size();
        double sum = c[n - 1];

        for (int i = n - 2; i >= 0; i--) {
            sum = sum * x + c[i];
        }

        return sum;
    }

    static double polyDiffEval(Vec c, const double x) {
        const int n = c.size();
        double sum = (n - 1) * c[n - 1];

        for (int i = n - 2; i >= 1; i--) {
            sum = sum * x + i * c[i];
        }

        return sum;
    }

    static CVec2 quadraticRoots(Vec c) {

        const double d = c[1] * c[1] - 4.0 * c[0] * c[2];

        if (d >= 0.0) {
            const double q = -0.5 * (c[1] + sgn(c[1]) * sqrt(d));

            return { q / c[2], c[0] / q };
        }
        
        const Complex q = -0.5 * (c[1] + sgn(c[1]) * sqrt(Complex(d)));

        return { q / c[2], c[0] / q };
        
    }

    static CVec3 cubicRoots(Vec c_) {
        constexpr double PI2 = 3.141592653589793238 * 2.0;
        //const double a_ = 1.0 / c_[3];
        const double a = c_[2] / c_[3], b = c_[1] / c_[3], c = c_[0] / c_[3];
        
        //constexpr double t = 1.0 / 3.0;
        const double a2 = a * a;
        const double Q = (a2 - 3.0 * b) / 9.0;
        const double R = (2.0 * a * a2 - 9.0 * a * b + 27.0 * c) / 54.0;
        const double R2 = R * R;
        const double Q3 = Q * Q * Q;
        const double k2 = -a / 3.0;
        
        if (R2 < Q3) {
            const double theta = acos(R / sqrt(Q3));
            const double k1 = -2.0 * sqrt(Q);

            return { k1 * cos(theta / 3.0) + k2, k1 * cos((theta + PI2) / 3.0) + k2, k1 * cos((theta - PI2) / 3.0) + k2 };
        }
        
        const double A = -sgn(R) * pow((fabs(R) + sqrt(R2 - Q3)), 1.0 / 3.0);

        if (fabs(A) < 1e-10) {
            return { k2, k2, k2 };
        }

        const double B = Q / A;
        const double p = A + B;
        const double m = (sqrt(3.0) / 2.0) * (A - B);

        return { p + k2, Complex(-0.5 * p + k2, m) ,  Complex(-0.5 * p + k2, -m)};
        
    }

    static Vec minimaxPoly(std::function<double(double)> func, Vec2 I, const int n, const double tol = 1e-7){
        
        double h = std::numeric_limits<double>::max();

        Mat A = constructMinimaxMat(I, n);
        Vec b, c, roots;
        Vec errs;

        b.resize(n + 2);
        roots.resize(n + 1);
        errs.resize(n);

        auto rootFunc = [&c, func](double x) { return func(x) - polyEval(c, x); };
        auto minFunc = [&c, func](double x) { return -fabs(func(x) - polyEval(c, x)); };

        int iter = 0;
        constexpr int MAX_ITER = 1000;

        while (iter < MAX_ITER) {
            iter++;

            // update matrix equation
            for (int i = 0; i < n + 2; i++) {
                double s = A[i][1];

                b[i] = (func(s));

                for (int j = 2; j < n + 1; j++) {
                    s *= A[i][1];
                    A[i][j] = s;
                }
            }

            // solve system
            c = linSolve(A, b);
            
            ////check if converged
            //if (fabs(h - c[n + 1]) < tol) {
            //    break;
            //}
            
            //remove h from coeffs list
            h = c[n + 1];
            c.pop_back();
            
            // find roots
            for (int i = 0; i < n + 1; i++) {
                roots[i] = findRoot(rootFunc, { A[i][1], A[i + 1][1] }, tol);
            }

            // find locations of worst case errors
            for (int i = 0; i < n; i++) {
                Vec2 m = minimize(minFunc, { roots[i], roots[i + 1] }, tol);

                A[i + 1][1] = m[0];
                errs[i] = fabs(m[1]);
            }

            double errMean = 0.0;

            for (const auto& err : errs) {
				errMean += err;
			}

            errMean /= n;

            double maxErr = 0.0;

            for (const auto& err : errs) {
                
                maxErr = std::max(maxErr, fabs(err - errMean));
            }

            c.push_back(errMean);

            if (maxErr < tol * errMean) {
				break;
			}

        }

        if (iter == MAX_ITER) {
            std::cout << "Failed to converge" << std::endl;
		}
         
        return c;
    }
};

