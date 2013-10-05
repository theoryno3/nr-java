
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.chapter09.Static.*;
import static com.snuggy.nr.chapter14.Static.*;
import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter06.*;
import com.snuggy.nr.chapter10.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class Fitexy {

    // Object for tting a straight line a C bx to a set of points .xi ; yi
    // / with errors in both xi and yi , respectively sigx and sigy. Call the
    // constructor to calculate the t. The answers are then available as the
    // variables a, b, siga, sigb, chi2, and q. Output quantities a and b
    // make y D aCbx minimize 2, whose value is returned as chi2. The 2
    // probability is returned as q, a small value indicating a poor t
    // (sometimes indicating underestimated errors). The standard errors on
    // a and b, siga and sigb, are not meaningful if either (i) the t is poor,
    // or (ii) b is so large that the data are consistent with a vertical
    // (in nite b) line. If siga and sigb are returned as BIG, then the data
    // are consistent with all values of b.
    
    private double a, b, siga, sigb, chi2, q; // Answers.
    private int ndat;
    private final double[] xx, yy, sx, sy, ww; // Variables that communicate with
                                         // Chixy.
    private $double aa_ref = $(0.0);
    private $double offs_ref = $(0.0);
    
    public double a() {
        return a;
    }
    
    public double b() {
        return b;
    }
    
    public double chi2() {
        return chi2;
    }
    
    public double q() {
        return q;
    }
    
    public double siga() {
        return siga;
    }
    
    public double sigb() {
        return sigb;
    }

    public Fitexy(final double[] x, final double[] y, final double[] sigx, final double[] sigy) throws NRException {
        // Constructor. Call with the input data x[0..ndat-1], y[0..ndat-1],
        // sigx[0..ndat-1], and sigy[0..ndat-1].
        ndat = (x.length);
        xx = doub_vec(ndat);
        yy = doub_vec(ndat);
        sx = doub_vec(ndat);
        sy = doub_vec(ndat);
        ww = doub_vec(ndat);
        final double POTN = 1.571000, BIG = 1.0e30, ACC = 1.0e-6;
        final double PI = 3.141592653589793238;
        Gamma gam = new Gamma();
        Brent brent = new Brent(ACC);
        Chixy chixy = new Chixy(xx, yy, sx, sy, ww, aa_ref, offs_ref); // Instantiate
                                                                       // a
                                                                       // Chixy
                                                                       // and
                                                                       // bind
                                                                       // it to
        int j; // our variables.
        double amx, amn, varx_ref[] = doub_ref(), vary_ref[] = doub_ref(), ang[] = new double[7], ch[] = new double[7], scale, bmn, bmx, d1_ref[] = doub_ref(), d2_ref[] = doub_ref(), r2, dum1_ref[] = doub_ref();
        avevar(x, dum1_ref, varx_ref); // Find the x and y variances, and scale
        // the data into the global variables for communication with the
        // function chixy.
        avevar(y, dum1_ref, vary_ref);
        scale = sqrt(varx_ref[0] / vary_ref[0]);
        for (j = 0; j < ndat; j++) {
            xx[j] = x[j];
            yy[j] = y[j] * scale;
            sx[j] = sigx[j];
            sy[j] = sigy[j] * scale;
            ww[j] = sqrt(SQR(sx[j]) + SQR(sy[j])); // Use both x and y weights
                                                   // in rst
        } // trial t.
        Fitab fit = new Fitab(xx, yy, ww);
        b = fit.b(); // Trial t for b.
        offs_ref.$(ang[0] = 0.0); // Construct several angles for reference
        // points, and make b an angle.
        ang[1] = atan(b);
        ang[3] = 0.0;
        ang[4] = ang[1];
        ang[5] = POTN;
        for (j = 3; j < 6; j++)
            ch[j] = chixy.eval(ang[j]);
        // Bracket the 2 minimum and then locate it with brent.
        brent.bracket(ang[0], ang[1], chixy);
        ang[0] = brent.ax();
        ang[1] = brent.bx();
        ang[2] = brent.cx();
        ch[0] = brent.fa();
        ch[1] = brent.fb();
        ch[2] = brent.fc();
        b = brent.minimize(chixy);
        chi2 = chixy.eval(b);
        a = aa_ref.$();
        q = gam.gammq(0.5 * (ndat - 2), chi2 * 0.5); // Compute 2 probability.
        r2 = 0.0;
        for (j = 0; j < ndat; j++)
            r2 += ww[j]; // Save the inverse sum of weights at
        r2 = 1.0 / r2; // the minimum.
        bmx = bmn = BIG; // Now, nd standard errors for b as
        // points where 2 offs=chi2+1.0; D 1.
        for (j = 0; j < 6; j++) { // Go through saved values to bracket
            // the desired roots. Note periodicity in slope angles.
            if (ch[j] > offs_ref.$()) {
                d1_ref[0] = abs(ang[j] - b);
                while (d1_ref[0] >= PI)
                    d1_ref[0] -= PI;
                d2_ref[0] = PI - d1_ref[0];
                if (ang[j] < b)
                    SWAP(d1_ref, d2_ref);
                if (d1_ref[0] < bmx)
                    bmx = d1_ref[0];
                if (d2_ref[0] < bmn)
                    bmn = d2_ref[0];
            }
        }
        if (bmx < BIG) { // Call zbrent to nd the roots.
            bmx = zbrent(chixy, b, b + bmx, ACC) - b;
            amx = aa_ref.$() - a;
            bmn = zbrent(chixy, b, b - bmn, ACC) - b;
            amn = aa_ref.$() - a;
            sigb = sqrt(0.5 * (bmx * bmx + bmn * bmn)) / (scale * SQR(cos(b)));
            siga = sqrt(0.5 * (amx * amx + amn * amn) + r2) / scale; // Error in
                                                                     // a has
                                                                     // additional
                                                                     // piece
        } else
            sigb = siga = BIG; // r2.
        a /= scale; // Unscale the answers.
        b = tan(b) / scale;
    }

    public void SWAP(double x_ref[], double y_ref[]) {
        double t = x_ref[0];
        x_ref[0] = y_ref[0];
        y_ref[0] = t;
    }

    class Chixy implements Func_Doub_To_Doub {
        // Captive functor of Fitexy, returns the value of .2 offs/ for the
        // slope b=tan(bang). Scaled data and offs are communicated via bound
        // references.
        private final double[] xx, yy, sx, sy, ww;
        private $double aa_ref, offs_ref;

        public Chixy(final double[] xxx, final double[] yyy, final double[] ssx, final double[] ssy,
                final double[] www, final $double aaa_ref, final $double ooffs_ref) {
            xx = (xxx);
            yy = (yyy);
            sx = (ssx);
            sy = (ssy);
            ww = (www);
            aa_ref = aaa_ref;
            offs_ref = (ooffs_ref);
        }

        // Constructor. Bind references back to Fitexy.

        public double eval(final double bang) {
            // The function as seen by Brent and zbrent.
            final double BIG = 1.0e30;
            int j, nn = xx.length;
            double ans, avex = 0.0, avey = 0.0, sumw = 0.0, b;
            b = tan(bang);
            for (j = 0; j < nn; j++) {
                ww[j] = SQR(b * sx[j]) + SQR(sy[j]);
                sumw += (ww[j] = (ww[j] < 1.0 / BIG ? BIG : 1.0 / ww[j]));
                avex += ww[j] * xx[j];
                avey += ww[j] * yy[j];
            }
            avex /= sumw;
            avey /= sumw;
            aa_ref.$(avey - b * avex);
            for (ans = -offs_ref.$(), j = 0; j < nn; j++)
                ans += ww[j] * SQR(yy[j] - aa_ref.$() - b * xx[j]);
            return ans;
        }
    }

}
