
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Complex.*;
import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;

import static java.lang.Math.*;

import java.lang.reflect.*;

import com.snuggy.nr.chapter17.*;
import com.snuggy.nr.util.*;

public class Static {

    public static double plegendre(final int l, final int m, final double x) throws NRException {
        // Computes the renormalized associated Legendre polynomial Pzm
        // l .x/, equation (6.7.8). Here m and l are integers satisfying
        // 0  m  l, while x lies in the range 1  x  1.
        final double PI = 3.141592653589793;
        int i, ll;
        double fact, oldfact, pll = 0.0, pmm, pmmp1, omx2;
        if (m < 0 || m > l || abs(x) > 1.0)
            throw new NRException("Bad arguments in routine plgndr");
        pmm = 1.0; // Compute Pzm m .
        if (m > 0) {
            omx2 = (1.0 - x) * (1.0 + x);
            fact = 1.0;
            for (i = 1; i <= m; i++) {
                pmm *= omx2 * fact / (fact + 1.0);
                fact += 2.0;
            }
        }
        pmm = sqrt((2 * m + 1) * pmm / (4.0 * PI));
        if ((m & 1) != 0)
            pmm = -pmm;
        if (l == m)
            return pmm;
        else { // Compute Pzm mC1.
            pmmp1 = x * sqrt(2.0 * m + 3.0) * pmm;
            if (l == (m + 1))
                return pmmp1;
            else { // Compute Pzm l , l >mC1.
                oldfact = sqrt(2.0 * m + 3.0);
                for (ll = m + 2; ll <= l; ll++) {
                    fact = sqrt((4.0 * ll * ll - 1.0) / (ll * ll - m * m));
                    pll = (x * pmmp1 - pmm / oldfact) * fact;
                    oldfact = fact;
                    pmm = pmmp1;
                    pmmp1 = pll;
                }
                return pll;
            }
        }
    }

    private static final double[] cof = { 57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
            -0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4, -.983744753048795646e-4,
            .158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
            .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5 };

    public static double gammln(final double xx) throws NRException {
        // Returns the value lnOE
        // .xx/ for xx > 0.
        int j;
        double x, tmp, y, ser;
        if (xx <= 0)
            throw new NRException("bad arg in gammln");
        y = x = xx;
        tmp = x + 5.24218750000000000; // Rational 671/128.
        tmp = (x + 0.5) * log(tmp) - tmp;
        ser = 0.999999999999997092;
        for (j = 0; j < 14; j++)
            ser += cof[j] / ++y;
        return tmp + log(2.5066282746310005 * ser / x);
    }

    private static final double[] factrl_a = doub_vec(171);
    private static boolean factrl_init = true;

    public static double factrl(final int n) throws NRException {
        // Returns the value nä as a floating-point number.
        if (factrl_init) {
            factrl_init = false;
            factrl_a[0] = 1.;
            for (int i = 1; i < 171; i++)
                factrl_a[i] = i * factrl_a[i - 1];
        }
        if (n < 0 || n > 170)
            throw new NRException("factrl out of range");
        return factrl_a[n];
    }

    private static final int NTOP = 2000;
    private static final double[] a = doub_vec(NTOP);
    private static boolean factln_init = true;

    public static double factln(final int n) throws NRException {
        // Returns ln.nä/.
        if (factln_init) {
            factln_init = false;
            for (int i = 0; i < NTOP; i++)
                a[i] = gammln(i + 1.);
        }
        if (n < 0)
            throw new NRException("negative arg in factln");
        if (n < NTOP)
            return a[n];
        return gammln(n + 1.); // Out of range of table.
    }
    
    public static double bico(final int n, final int k) throws NRException {
        // Returns the binomial coefficient
        // n
        // k as a floating-point number.
        if (n < 0 || k < 0 || k > n)
            throw new NRException("bad args in bico");
        if (n < 171)
            return floor(0.5 + factrl(n) / (factrl(k) * factrl(n - k)));
        return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
        // The floor function cleans up roundoff error for smaller values of n
        // and k.
    }
    
    public static double beta(final double z, final double w) throws NRException {
	    return exp(gammln(z)+gammln(w)-gammln(z+w));
	}
    
    private static final int expint_MAXIT = 100;
    private static final double expint_EULER = 0.577215664901533, expint_EPS = EPS(), // numeric_limits<Doub>::epsilon(),
            expint_BIG = Double.MAX_VALUE * expint_EPS; // numeric_limits<Doub>::max()*EPS;

    public static double expint(final int n, final double x) throws NRException {
        // Evaluates the exponential integral En.x/.
        // Here MAXIT is the maximum allowed number of iterations; EULER
        // is Eulerís constant ; EPS is the desired relative error, not
        // smaller than the machine precision; BIG is a number near the
        // largest representable floating-point number.
        int i, ii, nm1 = n - 1;
        double a, b, c, d, del, fact, h, psi, ans;
        if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
            throw new NRException("bad arguments in expint");
        if (n == 0)
            ans = exp(-x) / x; // Special case.
        else {
            if (x == 0.0)
                ans = 1.0 / nm1; // Another special case.
            else {
                if (x > 1.0) { // Lentzís algorithm (5.2).
                    b = x + n;
                    c = expint_BIG;
                    d = 1.0 / b;
                    h = d;
                    for (i = 1; i <= expint_MAXIT; i++) {
                        a = -i * (nm1 + i);
                        b += 2.0;
                        d = 1.0 / (a * d + b); // Denominators cannot be zero.
                        c = b + a / c;
                        del = c * d;
                        h *= del;
                        if (abs(del - 1.0) <= expint_EPS) {
                            ans = h * exp(-x);
                            return ans;
                        }
                    }
                    throw new NRException("continued fraction failed in expint");
                } else { // Evaluate series.
                    ans = (nm1 != 0 ? 1.0 / nm1 : -log(x) - expint_EULER); // Set
                                                                           // first
                                                                           // term.
                    fact = 1.0;
                    for (i = 1; i <= expint_MAXIT; i++) {
                        fact *= -x / i;
                        if (i != nm1)
                            del = -fact / (i - nm1);
                        else {
                            psi = -expint_EULER; // Compute .n/.
                            for (ii = 1; ii <= nm1; ii++)
                                psi += 1.0 / ii;
                            del = fact * (-log(x) + psi);
                        }
                        ans += del;
                        if (abs(del) < abs(ans) * expint_EPS)
                            return ans;
                    }
                    throw new NRException("series failed in expint");
                }
            }
        }
        return ans;
    }

    private static final int ei_MAXIT = 100;
    private static final double er_EULER = 0.577215664901533, ei_EPS = EPS(), // numeric_limits<Doub>::epsilon(),
            ei_FPMIN = Double.MIN_VALUE / ei_EPS; // numeric_limits<Doub>::min()/EPS;

    public static double ei(final double x) throws NRException {
        // Computes the exponential integral Ei.x/ for x > 0.
        // Here MAXIT is the maximum number of iterations allowed; EULER
        // is Eulerís constant ; EPS is the relative error, or absolute
        // error near the zero of Ei at x D 0:3725; FPMIN is a number close
        // to the smallest representable floating-point number.
        int k;
        double fact, prev, sum, term;
        if (x <= 0.0)
            throw new NRException("Bad argument in ei");
        if (x < ei_FPMIN)
            return log(x) + er_EULER; // Special case: Avoid failure of
                                      // convergence
        if (x <= -log(ei_EPS)) { // test because of underflow.
            sum = 0.0; // Use power series.
            fact = 1.0;
            for (k = 1; k <= ei_MAXIT; k++) {
                fact *= x / k;
                term = fact / k;
                sum += term;
                if (term < ei_EPS * sum)
                    break;
            }
            if (k > ei_MAXIT)
                throw new NRException("Series failed in ei");
            return sum + log(x) + er_EULER;
        } else { // Use asymptotic series.
            sum = 0.0; // Start with second term.
            term = 1.0;
            for (k = 1; k <= ei_MAXIT; k++) {
                prev = term;
                term *= k / x;
                if (term < ei_EPS)
                    break;
                // Since final sum is greater than one, term itself approximates
                // the relative error.
                if (term < prev)
                    sum += term; // Still converging: Add new term.
                else {
                    sum -= prev; // Diverging: Subtract previous term and
                    break; // exit.
                }
            }
            return exp(x) * (1.0 + sum) / x;
        }
    }

    public static double plgndr(final int l, final int m, final double x) throws NRException {
        // Computes the associated Legendre polynomial Pm
        // l .x/, equation (6.7.4). Here m and l are integers satisfying
        // 0  m  l, while x lies in the range 1  x  1. These functions will
        // overflow for m & 80.
        final double PI = 3.141592653589793238;
        if (m < 0 || m > l || abs(x) > 1.0)
            throw new NRException("Bad arguments in routine plgndr");
        double prod = 1.0;
        for (int j = l - m + 1; j <= l + m; j++)
            prod *= j;
        return sqrt(4.0 * PI * prod / (2 * l + 1)) * plegendre(l, m, x);
    }

    private static final int frenel_MAXIT = 100;
    private static final double frenel_PI = 3.141592653589793238, frenel_PIBY2 = (frenel_PI / 2.0), frenel_XMIN = 1.5,
            frenel_EPS = EPS(), // numeric_limits<Doub>::epsilon(),
            frenel_FPMIN = Double.MIN_VALUE, // numeric_limits<Doub>::min(),
            frenel_BIG = Double.MAX_VALUE * frenel_EPS; // numeric_limits<Doub>::max()*EPS;

    public static Complex frenel(final double x) throws NRException {
        // Computes the Fresnel integrals S.x/ and C.x/ for all real x. C.x/
        // is returned as the real part of cs and S.x/ as the imaginary part.
        // Here MAXIT is the maximum number of iterations allowed; EPS is the
        // relative error;
        // FPMIN is a number near the smallest representable floating-point
        // number; BIG is a
        // number near the machine overflow limit; and XMIN is the dividing line
        // between using
        // the series and continued fraction.
        boolean odd;
        int k, n;
        double a, ax, fact, pix2, sign, sum, sumc, sums, term, test;
        Complex b, cc, d, h, del, cs;
        if ((ax = abs(x)) < sqrt(frenel_FPMIN)) { // Special case: Avoid failure
                                                  // of
            // convergence
            cs = complex(ax); // test because of underflow.
        } else if (ax <= frenel_XMIN) { // Evaluate both series simultaneously.
            sum = sums = 0.0;
            sumc = ax;
            sign = 1.0;
            fact = frenel_PIBY2 * ax * ax;
            odd = true;
            term = ax;
            n = 3;
            for (k = 1; k <= frenel_MAXIT; k++) {
                term *= fact / k;
                sum += sign * term / n;
                test = abs(sum) * frenel_EPS;
                if (odd) {
                    sign = -sign;
                    sums = sum;
                    sum = sumc;
                } else {
                    sumc = sum;
                    sum = sums;
                }
                if (term < test)
                    break;
                odd = !odd;
                n += 2;
            }
            if (k > frenel_MAXIT)
                throw new NRException("series failed in frenel");
            cs = complex(sumc, sums);
        } else { // Evaluate continued fraction by modified
            pix2 = frenel_PI * ax * ax; // Lentzís method (5.2).
            b = complex(1.0, -pix2);
            cc = complex(frenel_BIG);
            // d=h=1.0/b;
            d = h = divide(1.0, b);
            n = -1;
            for (k = 2; k <= frenel_MAXIT; k++) {
                n += 2;
                a = -n * (n + 1);
                // b += 4.0;
                b = plus(b, 4.0);
                // d=1.0/(a*d+b); // Denominators cannot be zero.
                d = divide(1.0, plus(times(a, d), b));
                // cc=b+a/cc;
                cc = plus(b, divide(a, cc));
                // del=cc*d;
                del = times(cc, d);
                // h *= del;
                h = times(h, del);
                if (abs(real(del) - 1.0) + abs(imag(del)) <= frenel_EPS)
                    break;
            }
            if (k > frenel_MAXIT)
                throw new NRException("cf failed in frenel");
            // h *= Complex(ax,-ax);
            h = times(h, complex(ax, -ax));
            // cs=Complex(0.5,0.5)
            // *(1.0-Complex(cos(0.5*pix2),sin(0.5*pix2))*h);
            cs = times(complex(0.5, 0.5), minus(1.0, times(complex(cos(times(0.5, pix2)), sin(times(0.5, pix2))), h)));
        }
        if (x < 0.0)
            cs = minus(cs); // Use antisymmetry.
        return cs;
    }

    private static final int MAXIT = 100; // Maximum number of iterations
                                          // allowed.
    private static final double cisi_EULER = 0.577215664901533, cisi_PIBY2 = 1.570796326794897, cisi_TMIN = 2.0,
            cisi_EPS = EPS(), // numeric_limits<Doub>::epsilon(),
            cisi_FPMIN = Double.MIN_VALUE * 4.0, // numeric_limits<Doub>::min()*4.0,
            cisi_BIG = Double.MAX_VALUE * cisi_EPS; // numeric_limits<Doub>::max()*EPS;

    // Here EULER is Eulerís constant ; PIBY2 is =2; TMIN is the dividing
    // line between using the series and continued fraction; EPS is the
    // relative error, or absolute error near a zero of Ci.x/;
    // FPMIN is a number close to the smallest representable floating-point
    // number; and BIG is a number near the machine overflow limit.

    public static Complex cisi(final double x) throws NRException {
        // Computes the cosine and sine integrals Ci.x/ and Si.x/.
        // The function Ci.x/ is returned as the real part of cs, and Si.x/
        // as the imaginary part. Ci.0/ is returned as a large negative
        // number and no error message is generated. For x < 0 the routine
        // returns Ci.x/ and you must supply the i yourself.
        int i, k;
        boolean odd;
        double a, err, fact, sign, sum, sumc, sums, t, term;
        Complex h, b, c, d, del, cs;
        if ((t = abs(x)) == 0.0)
            return complex(-cisi_BIG); // Special case.
        if (t > cisi_TMIN) { // Evaluate continued fraction by modified
            b = complex(1.0, t); // Lentzís method (5.2).
            c = complex(cisi_BIG, 0.0);
            // d=h=1.0/b;
            d = h = divide(1.0, b);
            for (i = 1; i < MAXIT; i++) {
                a = -i * i;
                // b += 2.0;
                b = plus(b, 2.0);
                // d=1.0/(a*d+b); // Denominators cannot be zero.
                d = divide(1.0, plus(times(a, d), b));
                // c=b+a/c;
                c = plus(b, divide(a, c));
                // del=c*d;
                del = times(c, d);
                // h *= del;
                h = times(h, del);
                if (abs(real(del) - 1.0) + abs(imag(del)) <= cisi_EPS)
                    break;
            }
            if (i >= MAXIT)
                throw new NRException("cf failed in cisi");
            // h=Complex(cos(t),-sin(t))*h;
            h = times(complex(cos(t), -sin(t)), h);
            // cs= -conj(h)+Complex(0.0,PIBY2);
            cs = plus(minus(conj(h)), complex(0.0, cisi_PIBY2));
        } else { // Evaluate both series simultaneously.
            if (t < sqrt(cisi_FPMIN)) { // Special case: Avoid failure of
                                        // convergence
                sumc = 0.0; // test because of underflow.
                sums = t;
            } else {
                sum = sums = sumc = 0.0;
                sign = fact = 1.0;
                odd = true;
                for (k = 1; k <= MAXIT; k++) {
                    fact *= t / k;
                    term = fact / k;
                    sum += sign * term;
                    err = term / abs(sum);
                    if (odd) {
                        sign = -sign;
                        sums = sum;
                        sum = sumc;
                    } else {
                        sumc = sum;
                        sum = sums;
                    }
                    if (err < cisi_EPS)
                        break;
                    odd = !odd;
                }
                if (k > MAXIT)
                    throw new NRException("maxits exceeded in cisi");
            }
            cs = complex(sumc + log(t) + cisi_EULER, sums);
        }
        if (x < 0.0)
            cs = conj(cs);
        return cs;
    }

    private static final int NMAX = 6;
    private static final double[] c = doub_vec(NMAX);
    private static boolean init = true;
    private static final double H = 0.4, A1 = 2.0 / 3.0, A2 = 0.4, A3 = 2.0 / 7.0;

    public static double dawson(final double x) {
        // Returns Dawsonís integral F.x/ D exp.x2/
        // R x
        // 0 exp.t2/dt for any real x.
        int i, n0; // Flag is true if we need to initialize, else false.
        double d1, d2, e1, e2, sum, x2, xp, xx, ans;
        if (init) {
            init = false;
            for (i = 0; i < NMAX; i++)
                c[i] = exp(-SQR((2.0 * i + 1.0) * H));
        }
        if (abs(x) < 0.2) { // Use series expansion.
            x2 = x * x;
            ans = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
        } else { // Use sampling theorem representation.
            xx = abs(x);
            n0 = 2 * Int(0.5 * xx / H + 0.5);
            xp = xx - n0 * H;
            e1 = exp(2.0 * xp * H);
            e2 = e1 * e1;
            d1 = n0 + 1;
            d2 = d1 - 2.0;
            sum = 0.0;
            for (i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
                sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
            ans = 0.5641895835 * SIGN(exp(-xp * xp), x) * sum; // Constant is 1=
                                                               // p .
        }
        return ans;
    }

    public static double invxlogx(final double y) throws NRException {
        // For negative y, 0 > y > e1, return x such that y D x log.x/. The
        // value returned is always the smaller of the two roots and is in the
        // range 0 < x < e1.
        final double ooe = 0.367879441171442322;
        double t, u, to = 0.;
        if (y >= 0. || y <= -ooe)
            throw new NRException("no such inverse value");
        if (y < -0.2)
            u = log(ooe - sqrt(2 * ooe * (y + ooe))); // First approximation by
                                                      // inverse
        else
            u = -10.; // of Taylor series.
        do { // See text for derivation.
            u += (t = (log(y / u) - u) * (u / (1. + u)));
            if (t < 1.e-8 && abs(t + to) < 0.01 * abs(t))
                break;
            to = t;
        } while (abs(t / u) > 1.e-15);
        return exp(u);
    }

    private static final double rf_ERRTOL = 0.0025, rf_THIRD = 1.0 / 3.0, rf_C1 = 1.0 / 24.0, rf_C2 = 0.1,
            rf_C3 = 3.0 / 44.0, rf_C4 = 1.0 / 14.0;
    private static final double rf_TINY = 5.0 * Double.MIN_VALUE, // numeric_limits<Doub>::min(),
            rf_BIG = 0.2 * Double.MAX_VALUE; // numeric_limits<Doub>::max();

    public static double rf(final double x, final double y, final double z) throws NRException {
        // Computes Carlsonís elliptic integral of the first kind,
        // RF .x; y; z/. x, y, and z must be nonnegative, and at most one
        // can be zero.
        double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;
        if (MIN(MIN(x, y), z) < 0.0 || MIN(MIN(x + y, x + z), y + z) < rf_TINY || MAX(MAX(x, y), z) > rf_BIG)
            throw new NRException("invalid arguments in rf");
        xt = x;
        yt = y;
        zt = z;
        do {
            sqrtx = sqrt(xt);
            sqrty = sqrt(yt);
            sqrtz = sqrt(zt);
            alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
            xt = 0.25 * (xt + alamb);
            yt = 0.25 * (yt + alamb);
            zt = 0.25 * (zt + alamb);
            ave = rf_THIRD * (xt + yt + zt);
            delx = (ave - xt) / ave;
            dely = (ave - yt) / ave;
            delz = (ave - zt) / ave;
        } while (MAX(MAX(abs(delx), abs(dely)), abs(delz)) > rf_ERRTOL);
        e2 = delx * dely - delz * delz;
        e3 = delx * dely * delz;
        return (1.0 + (rf_C1 * e2 - rf_C2 - rf_C3 * e3) * e2 + rf_C4 * e3) / sqrt(ave);
    }

    private static final double rd_ERRTOL = 0.0015, rd_C1 = 3.0 / 14.0, rd_C2 = 1.0 / 6.0, rd_C3 = 9.0 / 22.0,
            rd_C4 = 3.0 / 26.0, rd_C5 = 0.25 * rd_C3, rd_C6 = 1.5 * rd_C4;
    private static final double rd_TINY = 2.0 * pow(Double.MAX_VALUE/*
                                                                     * numeric_limits
                                                                     * <
                                                                     * Doub>::max
                                                                     * ()
                                                                     */, -2. / 3.), rd_BIG = 0.1 * rd_ERRTOL
            * pow(Double.MIN_VALUE/* numeric_limits<Doub>::min() */, -2. / 3.);

    public static double rd(final double x, final double y, final double z) throws NRException {
        // Computes Carlsonís elliptic integral of the second kind, RD.x;
        // y; z/. x and y must be nonnegative, and at most one can be zero.
        // z must be positive.
        double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;
        if (MIN(x, y) < 0.0 || MIN(x + y, z) < rd_TINY || MAX(MAX(x, y), z) > rd_BIG)
            throw new NRException("invalid arguments in rd");
        xt = x;
        yt = y;
        zt = z;
        sum = 0.0;
        fac = 1.0;
        do {
            sqrtx = sqrt(xt);
            sqrty = sqrt(yt);
            sqrtz = sqrt(zt);
            alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
            sum += fac / (sqrtz * (zt + alamb));
            fac = 0.25 * fac;
            xt = 0.25 * (xt + alamb);
            yt = 0.25 * (yt + alamb);
            zt = 0.25 * (zt + alamb);
            ave = 0.2 * (xt + yt + 3.0 * zt);
            delx = (ave - xt) / ave;
            dely = (ave - yt) / ave;
            delz = (ave - zt) / ave;
        } while (MAX(MAX(abs(delx), abs(dely)), abs(delz)) > rd_ERRTOL);
        ea = delx * dely;
        eb = delz * delz;
        ec = ea - eb;
        ed = ea - 6.0 * eb;
        ee = ed + ec + ec;
        return 3.0
                * sum
                + fac
                * (1.0 + ed * (-rd_C1 + rd_C5 * ed - rd_C6 * delz * ee) + delz
                        * (rd_C2 * ee + delz * (-rd_C3 * ec + delz * rd_C4 * ea))) / (ave * sqrt(ave));
    }

    private static final double rj_ERRTOL = 0.0015, rj_C1 = 3.0 / 14.0, rj_C2 = 1.0 / 3.0, rj_C3 = 3.0 / 22.0,
            rj_C4 = 3.0 / 26.0, rj_C5 = 0.75 * rj_C3, rj_C6 = 1.5 * rj_C4, rj_C7 = 0.5 * rj_C2, rj_C8 = rj_C3 + rj_C3;
    private static final double rj_TINY = pow(5.0 * Double.MIN_VALUE
    /* numeric_limits<Doub>::min() */, 1. / 3.), rj_BIG = 0.3 * pow(0.2 * Double.MAX_VALUE/*
                                                                                           * numeric_limits
                                                                                           * <
                                                                                           * Doub
                                                                                           * >
                                                                                           * :
                                                                                           * :
                                                                                           * max
                                                                                           * (
                                                                                           * )
                                                                                           */, 1. / 3.);

    public static double rj(final double x, final double y, final double z, final double p) throws NRException {
        // Computes Carlsonís elliptic integral of the third kind,
        // RJ .x; y; z;p/. x, y, and z must be nonnegative, and at most one
        // can be zero. p must be nonzero. If p < 0, the Cauchy principal value
        // is returned.

        double a = 0.0, alamb, alpha, ans, ave, b = 0.0, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt, rcx = 0.0, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt, zt;
        if (MIN(MIN(x, y), z) < 0.0 || MIN(MIN(x + y, x + z), MIN(y + z, abs(p))) < rj_TINY
                || MAX(MAX(x, y), MAX(z, abs(p))) > rj_BIG)
            throw new NRException("invalid arguments in rj");
        sum = 0.0;
        fac = 1.0;
        if (p > 0.0) {
            xt = x;
            yt = y;
            zt = z;
            pt = p;
        } else {
            xt = MIN(MIN(x, y), z);
            zt = MAX(MAX(x, y), z);
            yt = x + y + z - xt - zt;
            a = 1.0 / (yt - p);
            b = a * (zt - yt) * (yt - xt);
            pt = yt + b;
            rho = xt * zt / yt;
            tau = p * pt / yt;
            rcx = rc(rho, tau);
        }
        do {
            sqrtx = sqrt(xt);
            sqrty = sqrt(yt);
            sqrtz = sqrt(zt);
            alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
            alpha = SQR(pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
            beta = pt * SQR(pt + alamb);
            sum += fac * rc(alpha, beta);
            fac = 0.25 * fac;
            xt = 0.25 * (xt + alamb);
            yt = 0.25 * (yt + alamb);
            zt = 0.25 * (zt + alamb);
            pt = 0.25 * (pt + alamb);
            ave = 0.2 * (xt + yt + zt + pt + pt);
            delx = (ave - xt) / ave;
            dely = (ave - yt) / ave;
            delz = (ave - zt) / ave;
            delp = (ave - pt) / ave;
        } while (MAX(MAX(abs(delx), abs(dely)), MAX(abs(delz), abs(delp))) > rj_ERRTOL);
        ea = delx * (dely + delz) + dely * delz;
        eb = delx * dely * delz;
        ec = delp * delp;
        ed = ea - 3.0 * ec;
        ee = eb + 2.0 * delp * (ea - ec);
        ans = 3.0
                * sum
                + fac
                * (1.0 + ed * (-rj_C1 + rj_C5 * ed - rj_C6 * ee) + eb * (rj_C7 + delp * (-rj_C8 + delp * rj_C4)) + delp
                        * ea * (rj_C2 - delp * rj_C3) - rj_C2 * delp * ec) / (ave * sqrt(ave));
        if (p <= 0.0)
            ans = a * (b * ans + 3.0 * (rcx - rf(xt, yt, zt)));
        return ans;
    }

    private static final double rc_ERRTOL = 0.0012, rc_THIRD = 1.0 / 3.0, rc_C1 = 0.3, rc_C2 = 1.0 / 7.0,
            rc_C3 = 0.375, rc_C4 = 9.0 / 22.0;
    private static final double rc_TINY = 5.0 * Double.MIN_VALUE, // numeric_limits<Doub>::min(),
            rc_BIG = 0.2 * Double.MAX_VALUE, // numeric_limits<Doub>::max(),
            rc_COMP1 = 2.236 / sqrt(rc_TINY), rc_COMP2 = SQR(rc_TINY * rc_BIG) / 25.0;

    public static double rc(final double x, final double y) throws NRException {
        // Computes Carlsonís degenerate elliptic integral, RC .x; y/. x
        // must be nonnegative and y must be nonzero. If y < 0, the Cauchy
        // principal value is returned.
        double alamb, ave, s, w, xt, yt;
        if (x < 0.0 || y == 0.0 || (x + abs(y)) < rc_TINY || (x + abs(y)) > rc_BIG
                || (y < -rc_COMP1 && x > 0.0 && x < rc_COMP2))
            throw new NRException("invalid arguments in rc");
        if (y > 0.0) {
            xt = x;
            yt = y;
            w = 1.0;
        } else {
            xt = x - y;
            yt = -y;
            w = sqrt(x) / sqrt(xt);
        }
        do {
            alamb = 2.0 * sqrt(xt) * sqrt(yt) + yt;
            xt = 0.25 * (xt + alamb);
            yt = 0.25 * (yt + alamb);
            ave = rc_THIRD * (xt + yt + yt);
            s = (yt - ave) / ave;
        } while (abs(s) > rc_ERRTOL);
        return w * (1.0 + s * s * (rc_C1 + s * (rc_C2 + s * (rc_C3 + s * rc_C4)))) / sqrt(ave);
    }

    public static double ellf(final double phi, final double ak) throws NRException {
        // Legendre elliptic integral of the first kind F.;k/, evaluated
        // using Carlsonís function RF. The argument ranges are 0    =2, 0 
        // k sin   1.
        double s = sin(phi);
        return s * rf(SQR(cos(phi)), (1.0 - s * ak) * (1.0 + s * ak), 1.0);
    }

    public static double elle(final double phi, final double ak) throws NRException {
        // Legendre elliptic integral of the second kind E.;k/, evaluated
        // using Carlsonís functions RD and RF . The argument ranges
        // are 0    =2, 0  k sin   1.
        double cc, q, s;
        s = sin(phi);
        cc = SQR(cos(phi));
        q = (1.0 - s * ak) * (1.0 + s * ak);
        return s * (rf(cc, q, 1.0) - (SQR(s * ak)) * rd(cc, q, 1.0) / 3.0);
    }

    public static double ellpi(final double phi, final double en, final double ak) throws NRException {
        // Legendre elliptic integral of the third kind Ö.; n; k/, evaluated
        // using Carlsonís functions RJ and RF . (Note that the sign convention
        // on n is opposite that of Abramowitz and Stegun.)
        // The ranges of  and k are 0    =2, 0  k sin   1.
        double cc, enss, q, s;
        s = sin(phi);
        enss = en * s * s;
        cc = SQR(cos(phi));
        q = (1.0 - s * ak) * (1.0 + s * ak);
        return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
    }

    private static final double CA = 1.0e-8; // The accuracy is the square of
                                             // CA.

    public static void sncndn(final double uu, final double emmc, 
            final $double sn, final $double cn, final $double dn) {
        // Returns the Jacobian elliptic functions sn.u; kc/, cn.u; kc/,
        // and dn.u; kc/. Here uu D u, while emmc D k2 c .
        boolean bo;
        int i, ii, l = 0;
        double a, b, c = 0.0, d = 0.0, emc, u;
        final double[] em = doub_vec(13), en = doub_vec(13);
        emc = emmc;
        u = uu;
        if (emc != 0.0) {
            bo = (emc < 0.0);
            if (bo) {
                d = 1.0 - emc;
                emc /= -1.0 / d;
                u *= (d = sqrt(d));
            }
            a = 1.0;
            dn.$(1.0);
            for (i = 0; i < 13; i++) {
                l = i;
                em[i] = a;
                en[i] = (emc = sqrt(emc));
                c = 0.5 * (a + emc);
                if (abs(a - emc) <= CA * a)
                    break;
                emc *= a;
                a = c;
            }
            u *= c;
            sn.$(sin(u));
            cn.$(cos(u));
            if (sn.$() != 0.0) {
                a = cn.$() / sn.$();
                c *= a;
                for (ii = l; ii >= 0; ii--) {
                    b = em[ii];
                    a *= c;
                    c *= dn.$();
                    dn.$((en[ii] + a) / (b + a));
                    a = c / b;
                }
                a = 1.0 / sqrt(c * c + 1.0);
                sn.$((sn.$() >= 0.0 ? a : -a));
                cn.$(c * sn.$());
            }
            if (bo) {
                a = dn.$();
                dn.$(cn.$());
                cn.$(a);
                sn.$(sn.$() / d);
            }
        } else {
            cn.$(1.0 / cosh(u));
            dn.$(cn.$());
            sn.$(tanh(u));
        }
    }

    public static $$<Complex> hypgeo(final Complex a, final Complex b, final Complex c,
                                 final Complex z) throws NRException, InstantiationException, IllegalAccessException, NoSuchMethodException, SecurityException, IllegalArgumentException, InvocationTargetException {
        // Complex hypergeometric function 2F1 for complex a; b; c, and z, by
        // direct integration of the hypergeometric equation in the complex
        // plane. The branch cut is taken to lie along the real axis, Re z > 1.
        final double atol = 1.0e-14, rtol = 1.0e-14; // Accuracy parameters.
        $$<Complex> ans = $$(complex());
        Complex dz = complex(), z0 = complex();
        @SuppressWarnings("unchecked")
        $$<Complex>[] y = obj_vec_nulls($$.class, 2);
        y[0] = $$(complex());
        y[1] = $$(complex());
        final double[] yy = doub_vec(4);
        if (norm2(z) <= 0.25) { // Use series...
            hypser(a, b, c, z, ans, y[1]);
            return ans;
        }
        // ...or pick a starting point for the path integration.
        else if (real(z) < 0.0)
            z0 = complex(-0.5, 0.0);
        else if (real(z) <= 1.0)
            z0 = complex(0.5, 0.0);
        else
            z0 = complex(0.0, imag(z) >= 0.0 ? 0.5 : -0.5);
        // dz=z-z0;
        dz = minus(z, z0);
        hypser(a, b, c, z0, y[0], y[1]); // Get starting function and
                                         // derivative.
        yy[0] = real(y[0].$());
        yy[1] = imag(y[0].$());
        yy[2] = real(y[1].$());
        yy[3] = imag(y[1].$());
        Hypderiv d = new Hypderiv(a, b, c, z0, dz); // Set up the functor for
                                                    // the derivatives.
        
        Output out = new Output(); // Suppress output in Odeint.
        Odeint ode = new Odeint(StepperBS.class, yy, 0.0, 1.0, atol, rtol, 0.1, 0.0, out, d);

        // The arguments to Odeint are the vector of independent variables,
        // the starting and ending values of the dependent variable, the
        // accuracy parameters, an initial guess for the stepsize, a minimum
        // stepsize, and the names of the output object and the derivative
        // object. The integration is performed by the Bulirsch-Stoer
        // stepping routine.
        ode.integrate();
        $$(y[0], complex(yy[0], yy[1]));
        return y[0];
    }

    public static void hypser(final Complex a, final Complex b, final Complex c, final Complex z,
            final $$<Complex> series, final $$<Complex> deriv) throws NRException {
        // Returns the hypergeometric series 2F1 and its derivative, iterating
        // to machine accuracy. For jzj  1=2 convergence is quite rapid.
        Complex fac = complex(1.0);
        Complex temp = complex(fac);
        Complex aa = complex(a);
        Complex bb = complex(b);
        Complex cc = complex(c);
        for (int n = 1; n <= 1000; n++) {
            // fac *= ((aa * bb) / cc);
            fac = (times(fac, divide(times(aa, bb), cc)));
            // deriv += fac;
            $$(deriv, plus(deriv.$(), fac));
            // fac *= ((1.0 / n) * z);
            fac = (times(fac, times(divide(1.0, n), z)));
            // series = temp + fac;
            $$(series, plus(temp, fac));
            if (equal(series.$(), temp))
                return;
            // temp = series;
            temp = complex(series.$());
            // aa += 1.0;
            aa = plus(aa, 1.0);
            // bb += 1.0;
            bb = plus(bb, 1.0);
            // cc += 1.0;
            cc = plus(cc, 1.0);
        }
        throw new NRException("convergence failure in hypser");
    }

    // A lower-order Chebyshev approximation produces a very concise routine,
    // though with only about single precision accuracy:

    public static double erfcc(final double x) {
        // Returns the complementary error function erfc.x/ with fractional
        // error everywhere less than
        // 1:2  107.
        double t, z = abs(x), ans;
        t = 2. / (2. + z);
        ans = t
                * exp(-z
                        * z
                        - 1.26551223
                        + t
                        * (1.00002368 + t
                                * (0.37409196 + t
                                        * (0.09678418 + t
                                                * (-0.18628806 + t
                                                        * (0.27886807 + t
                                                                * (-1.13520398 + t
                                                                        * (1.48851587 + t
                                                                                * (-0.82215223 + t * 0.17087277)))))))));
        return (x >= 0.0 ? ans : 2.0 - ans);
    }

}
