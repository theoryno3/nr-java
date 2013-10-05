
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Curve_interp {

    // Object for interpolating a curve specified by n points in dim dimensions.
    private int dim, n, in;
    private boolean cls; // Set if a closed curve.
    private final double[][] pts;
    private final double[] s;
    private final double[] ans;
    private Spline_interp[] srp;

    public Curve_interp(final double[][] ptsin) throws NRException, InstantiationException, IllegalAccessException {
        this(ptsin, false);
    }

    public Curve_interp(final double[][] ptsin, final boolean close) 
            throws NRException, InstantiationException, IllegalAccessException {
        // Constructor. The n  dim matrix ptsin inputs the data points. Input
        // close as 0 for an open curve, 1 for a closed curve. (For a closed
        // curve, the last data point should not duplicate the first — the
        // algorithm will connect them.)

        n = (nrows(ptsin));
        dim = (ncols(ptsin));
        in = (close ? 2 * n : n);
        cls = (close);
        pts = doub_mat(dim, in);
        s = doub_arr(in);
        ans = doub_arr(dim);
        srp = obj_arr(Spline_interp.class, dim);

        int i, ii, im, j, ofs;
        double ss, soff, db, de;
        ofs = close ? n / 2 : 0; // The trick for closed curves 
                                 // is to duplicate half a
        // period at the beginning and end, and then
        // use the middle half of the resulting spline.
        s[0] = 0.;
        for (i = 0; i < in; i++) {
            ii = (i - ofs + n) % n;
            im = (ii - 1 + n) % n;
            for (j = 0; j < dim; j++)
                pts[j][i] = ptsin[ii][j]; // Store transpose.
            if (i > 0) { // Accumulate arc length.
                s[i] = s[i - 1] + rad(ptsin[ii], 0, ptsin[im], 0);
                if (s[i] == s[i - 1])
                    throw new NRException("error in Curve_interp");
                // Consecutive points may not be identical. For a closed curve,
                // the last data point should not duplicate the first.
            }
        }
        ss = close ? s[ofs + n] - s[ofs] : s[n - 1] - s[0]; // Rescale parameter
                                                            // so that the
        soff = s[ofs]; // interval [0,1] is the whole curve (or one period).
        for (i = 0; i < in; i++)
            s[i] = (s[i] - soff) / ss;
        for (j = 0; j < dim; j++) { // Construct the splines using endpoint
                                    // derivatives.
            db = in < 4 ? 1.e99 : fprime(s, 0, pts[j], 0, 1);
            de = in < 4 ? 1.e99 : fprime(s, in - 1, pts[j], in - 1, -1);
            srp[j] = new Spline_interp(s, pts[j], db, de);
        }
    }

    // ~Curve_interp() {for (Int j=0;j<dim;j++) delete srp[j];}

    public final double[] interp(double t) throws NRException {
        // Interpolate a point on the stored curve. The point is parameterized
        // by t, in the range [0,1]. For open curves, values of t outside this
        // range will return extrapolations (dangerous!). For closed curves,
        // t is periodic with period 1.
        if (cls)
            t = t - floor(t);
        for (int j = 0; j < dim; j++)
            ans[j] = srp[j].interp(t);
        return ans;
    }

    public double fprime(final double[] x_arr, final int x_off, final double[] y_arr, final int y_off, final int pm) {
        // Utility for estimating the derivatives at the endpoints. x and y
        // point to the abscissa and ordinate of the endpoint. If pm is C1,
        // points to the right will be used (left endpoint); if it is 1,
        // points to the left will be used (right endpoint). See text, below.
        double s1 = x_arr[x_off + 0] - x_arr[x_off + pm * 1], 
               s2 = x_arr[x_off + 0] - x_arr[x_off + pm * 2], 
               s3 = x_arr[x_off + 0] - x_arr[x_off + pm * 3], 
               s12 = s1 - s2, s13 = s1 - s3, s23 = s2 - s3;
        return -(s1 * s2 / (s13 * s23 * s3)) * y_arr[y_off + pm * 3] + 
                (s1 * s3 / (s12 * s2 * s23)) * y_arr[y_off + pm * 2] - 
                (s2 * s3 / (s1 * s12 * s13)) * y_arr[y_off + pm * 1] + 
                (1. / s1 + 1. / s2 + 1. / s3) * y_arr[y_off + 0];
    }

    public double rad(final double[] p1_arr, final int p1_off, 
                        final double[] p2_arr, final int p2_off) {
        // Euclidean distance.
        double sum = 0.0;
        for (int i = 0; i < dim; i++)
            sum += SQR(p1_arr[p1_off + i] - p2_arr[p2_off + i]);
        return sqrt(sum);
    }

}
