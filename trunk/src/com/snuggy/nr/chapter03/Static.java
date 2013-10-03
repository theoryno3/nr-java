
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Static {

    public static void polcoe(final double[] x, final double[] y, final double[] cof) {
        // Given arrays x[0..n-1] and y[0..n-1] containing a tabulated
        // function yi D f.xi /, this routine returns an array of coefficients
        // cof[0..n-1], such that yi D
        // Pn1
        // jD0 cofj x
        // j
        // i .

        int k, j, i, n = x.length;
        double phi, ff, b;
        double[] s = doub_arr(n);
        for (i = 0; i < n; i++)
            s[i] = cof[i] = 0.0;
        s[n - 1] = -x[0];
        for (i = 1; i < n; i++) { // Coefficients si of the master polynomial
                                  // P.x/ are
            for (j = n - 1 - i; j < n - 1; j++)
                // found by recurrence.
                s[j] -= x[i] * s[j + 1];
            s[n - 1] -= x[i];
        }
        for (j = 0; j < n; j++) {
            phi = n;
            for (k = n - 1; k > 0; k--)
                // The quantity phi D
                // Q
                // j¤k.xj  xk/ is found as a
                phi = k * s[k] + x[j] * phi; // derivative of P.xj /.
            ff = y[j] / phi;
            b = 1.0; // Coefficients of polynomials in each term of the
            // Lagrange formula are found by synthetic division of
            // P.x/ by .x  xj /. The solution ck is accumulated.
            for (k = n - 1; k >= 0; k--) {
                cof[k] += b * ff;
                b = s[k] + x[j] * b;
            }
        }
    }

    public static void polcof(final double[] xa, final double[] ya, final double[] cof) 
            throws NRException {
        // Given arrays xa[0..n-1] and ya[0..n-1] containing a tabulated
        // function yai D f.xai /, this routine returns an array of
        // coefficients cof[0..n-1], such that yai D
        // Pn1
        // jD0 cofj xa
        // j
        // i .

        int k, j, i, n = xa.length;
        double xmin;
        double[] x = doub_arr(n), y = doub_arr(n);
        for (j = 0; j < n; j++) {
            x[j] = xa[j];
            y[j] = ya[j];
        }
        for (j = 0; j < n; j++) { // Fill a temporary vector whose size
            // decreases as each coefficient is
            // found.
            double[] x_t = doub_arr(n - j), y_t = doub_arr(n - j);
            for (k = 0; k < n - j; k++) {
                x_t[k] = x[k];
                y_t[k] = y[k];
            }
            Poly_interp interp = new Poly_interp(x, y, n - j);
            cof[j] = interp.rawinterp(0, 0.); // Extrapolate to x D 0.
            xmin = 1.0e99;
            k = -1;
            for (i = 0; i < n - j; i++) { // Find the remaining xi of smallest
                if (abs(x[i]) < xmin) { // absolute value
                    xmin = abs(x[i]);
                    k = i;
                }
                if (x[i] != 0.0) // (meanwhile reducing all the terms)
                    y[i] = (y[i] - cof[j]) / x[i];
            }
            for (i = k + 1; i < n - j; i++) { // and eliminate it.
                y[i - 1] = y[i];
                x[i - 1] = x[i];
            }
        }
    }

    private static int wt_d[] = { 
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2,  0,  0,  -1, 
        0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 
        0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, -3, 3, 0, 0,
        -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, 9, -9, 9, -9, 6,
        3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2, -6, 6, -6, 6, -4, 
        -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2, 2, -2, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, -6, 6, -6, 6, -3, -3,
        3, 3, -4, 4, 2, -2, -2, -2, -1, -1, 4, -4, 4, -4, 
        2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1 
    };
    
    static int[][] wt = int_mat(16, 16, wt_d);

    @Deprecated @Broken
    public static void bcucof(final double[] y, final double[] y1, final double[] y2, final double[] y12, 
                            final double d1, final double d2, final double[][] c) {
        // Given arrays y[0..3], y1[0..3], y2[0..3], and y12[0..3], containing
        // the function, gradients, and cross-derivative at the four grid
        // points of a rectangular grid cell (numbered counterclockwise from
        // the lower left), and given d1 and d2, the length of the grid cell
        // in the 1 and 2 directions, this routine returns the table
        // c[0..3][0..3] that is used by routine bcuint for bicubic
        // interpolation.

        int l, k, j, i;
        double xx, d1d2 = d1 * d2;
        double[] cl = doub_arr(16), x = doub_arr(16);
        for (i = 0; i < 4; i++) { // Pack a temporary vector x.
            x[i] = y[i];
            x[i + 4] = y1[i] * d1;
            x[i + 8] = y2[i] * d2;
            x[i + 12] = y12[i] * d1d2;
        }
        for (i = 0; i < 16; i++) { // Matrix-multiply by the stored table.
            xx = 0.0;
            for (k = 0; k < 16; k++)
                xx += wt[i][k] * x[k];
            cl[i] = xx;
        }
        l = 0;
        for (i = 0; i < 4; i++)
            // Unpack the result into the output table.
            for (j = 0; j < 4; j++)
                c[i][j] = cl[l++];
    }

    // The implementation of equation (3.6.6), which performs a bicubic
    // interpolation, gives back the interpolated function value and the two
    // gradient values, and uses the above routine bcucof, is simply:

    @Deprecated @Broken
    public static void bcuint(final double[] y, final double[] y1, final double[] y2, final double[] y12, 
                                final double x1l, final double x1u, 
                                final double x2l, final double x2u, 
                                final double x1, final double x2, 
                                final double ansy_ref[], final double ansy1_ref[], 
                                final double ansy2_ref[]) 
            throws NRException {
        // Bicubic interpolation within a grid square. Input quantities
        // are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower
        // and upper coordinates of the grid square in the 1 direction;
        // x2l and x2u likewise for the 2 direction; and x1,x2, the coordinates
        // of the desired point for the interpolation. The interpolated
        // function value is returned as ansy, and the interpolated gradient
        // values as ansy1 and ansy2. This routine calls bcucof.
        int i;
        double t, u, d1 = x1u - x1l, d2 = x2u - x2l;
        double[][] c = doub_mat(4, 4);
        bcucof(y, y1, y2, y12, d1, d2, c); // Get the c’s.
        if (x1u == x1l || x2u == x2l)
            throw new NRException("Bad input in routine bcuint");
        t = (x1 - x1l) / d1; // Equation (3.6.4).
        u = (x2 - x2l) / d2;
        ansy_ref[0] = ansy2_ref[0] = ansy1_ref[0] = 0.0;
        for (i = 3; i >= 0; i--) { // Equation (3.6.6).
            ansy_ref[0] = t * ansy_ref[0] + ((c[i][3] * u + c[i][2]) * u + c[i][1]) * u + c[i][0];
            ansy2_ref[0] = t * ansy2_ref[0] + (3.0 * c[i][3] * u + 2.0 * c[i][2]) * u + c[i][1];
            ansy1_ref[0] = u * ansy1_ref[0] + (3.0 * c[3][i] * t + 2.0 * c[2][i]) * t + c[1][i];
        }
        ansy1_ref[0] /= d1;
        ansy2_ref[0] /= d2;
    }

}
