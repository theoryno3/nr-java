
package com.snuggy.nr.chapter14;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.chapter08.Static.*;
import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.refs.Refs.*;
import com.snuggy.nr.refs.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.chapter06.*;
import com.snuggy.nr.util.*;

public class Static {

    public static void moment(final double[] data, final $double ave, final $double adev,
            final $double sdev, final $double var, final $double skew, final $double curt)
            throws NRException {
        // Given an array of data[0..n-1], this routine returns its mean ave,
        // average deviation adev, standard deviation sdev, variance var,
        // skewness skew, and kurtosis curt.
        int j, n = data.length;
        double ep = 0.0, s, p;
        if (n <= 1)
            throw new NRException("n must be at least 2 in moment");
        s = 0.0; // First pass to get the mean.
        for (j = 0; j < n; j++)
            s += data[j];
        ave.$(s / n);
        $(adev, $(var, $(skew, $(curt, 0.0)))); // Second
                                                // pass to
                                                // get the
                                                // first
                                                // (absolute),
                                                // second,
        // third, and fourth moments of the deviation from the mean.
        for (j = 0; j < n; j++) {
            adev.$(adev.$() + abs(s = data[j] - ave.$()));
            ep += s;
            var.$(var.$() + (p = s * s));
            skew.$(skew.$() + (p *= s));
            curt.$(curt.$() + (p *= s));
        }
        adev.$(adev.$() / n);
        var.$((var.$() - ep * ep / n) / (n - 1)); // Corrected two-pass
                                                           // formula.
        sdev.$(sqrt(var.$())); // Put the pieces together according to
                                        // the con
        if (var.$() != 0.0) { // ventional definitions.
            skew.$(skew.$() / (n * var.$() * sdev.$()));
            curt.$(curt.$() / (n * var.$() * var.$()) - 3.0);
        } else
            throw new NRException("No skew/kurtosis when variance = 0 (in moment)");
    }

    public static void ttest(final double[] data1, final double[] data2, final $double t, final $double prob)
            throws NRException {
        // Given the arrays data1[0..n1-1] and data2[0..n2-1], returns
        // Student’s t as t, and its pvalue as prob, small values of prob
        // indicating that the arrays have significantly different means. The
        // data arrays are assumed to be drawn from populations with the same
        // true variance.
        Beta beta = new Beta();
        final $double var1 = $(0.0), var2 = $(0.0), ave1 = $(0.0), ave2 = $(0.0);
        double svar, df;
        int n1 = data1.length, n2 = data2.length;
        avevar(data1, ave1, var1);
        avevar(data2, ave2, var2);
        df = n1 + n2 - 2; // Degrees of freedom.
        svar = ((n1 - 1) * var1.$() + (n2 - 1) * var2.$()) / df; // Pooled
                                                                       // variance.
        t.$((ave1.$() - ave2.$()) / sqrt(svar * (1.0 / n1 + 1.0 / n2)));
        prob.$(beta.betai(0.5 * df, 0.5, df / (df + t.$() * t.$()))); // See
                                                                                  // equation
                                                                                  // (6.14.11).
    }

    public static void avevar(final double[] data, final $double ave, final $double var) {
        // Given array data[0..n-1], returns its mean as ave and its variance as
        // var.
        double s, ep;
        int j, n = data.length;
        ave.$(0.0);
        for (j = 0; j < n; j++)
            ave.$(ave.$() + data[j]);
        ave.$(ave.$() / n);
        var.$(ep = 0.0);
        for (j = 0; j < n; j++) {
            s = data[j] - ave.$();
            ep += s;
            var.$(var.$() + (s * s));
        }
        var.$((var.$() - ep * ep / n) / (n - 1)); // Corrected two-pass
                                                           // formula (14.1.8).
    }

    public static void tutest(final double[] data1, final double[] data2, final $double t, final $double prob)
            throws NRException {
        // Given the arrays data1[0..n1-1] and data2[0..n2-1], this routine
        // returns Student’s t as t, and its p-value as prob, small values of
        // prob indicating that the arrays have significantly different means.
        // The data arrays are allowed to be drawn from populations with
        // unequal variances.
        Beta beta = new Beta();
        final $double var1 = $(0.0), var2 = $(0.0), ave1 = $(0.0), ave2 = $(0.0);
        double df;
        int n1 = data1.length, n2 = data2.length;
        avevar(data1, ave1, var1);
        avevar(data2, ave2, var2);
        t.$((ave1.$() - ave2.$()) / sqrt(var1.$() / n1 + var2.$() / n2));
        df = SQR(var1.$() / n1 + var2.$() / n2)
                / (SQR(var1.$() / n1) / (n1 - 1) + SQR(var2.$() / n2) / (n2 - 1));
        prob.$(beta.betai(0.5 * df, 0.5, df / (df + SQR(t.$()))));
    }

    public static void tptest(final double[] data1, final double[] data2, final $double t, final $double prob)
            throws NRException {
        // Given the paired arrays data1[0..n-1] and data2[0..n-1], this routine
        // returns Student’s t for paired data as t, and its p-value as prob,
        // small values of prob indicating a significant difference of means.
        Beta beta = new Beta();
        int j, n = data1.length;
        double sd, df, cov = 0.0;
        final $double var1 = $(0.0), var2 = $(0.0), ave1 = $(0.0), ave2 = $(0.0);
        avevar(data1, ave1, var1);
        avevar(data2, ave2, var2);
        for (j = 0; j < n; j++)
            cov += (data1[j] - ave1.$()) * (data2[j] - ave2.$());
        cov /= (df = n - 1);
        sd = sqrt((var1.$() + var2.$() - 2.0 * cov) / n);
        t.$((ave1.$() - ave2.$()) / sd);
        prob.$(beta.betai(0.5 * df, 0.5, df / (df + t.$() * t.$())));
    }

    public static void ftest(final double[] data1, final double[] data2, final $double f, final $double prob)
            throws NRException {
        // Given the arrays data1[0..n1-1] and data2[0..n2-1], this routine
        // returns the value of f, and its p-value as prob. Small values of
        // prob indicate that the two arrays have significantly different
        // variances.
        Beta beta = new Beta();
        double df1, df2;
        final $double var1 = $(0.0), var2 = $(0.0), ave1 = $(0.0), ave2 = $(0.0);
        int n1 = data1.length, n2 = data2.length;
        avevar(data1, ave1, var1);
        avevar(data2, ave2, var2);
        if (var1.$() > var2.$()) { // Make F the ratio of the larger
                                         // variance to the smaller
            f.$(var1.$() / var2.$()); // one.
            df1 = n1 - 1;
            df2 = n2 - 1;
        } else {
            f.$(var2.$() / var1.$());
            df1 = n2 - 1;
            df2 = n1 - 1;
        }
        prob.$(2.0 * beta.betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f.$())));
        if (prob.$() > 1.0)
            prob.$(2. - prob.$());
    }

    public static void chsone(final double[] bins, final double[] ebins, final $double df,
            final $double chsq, final $double prob) throws NRException {
        chsone(bins, ebins, df, chsq, prob, 1);
    }

    public static void chsone(final double[] bins, final double[] ebins, final $double df,
            final $double chsq, final $double prob, final int knstrn) throws NRException {
        // Given the array bins[0..nbins-1] containing the observed numbers
        // of events, and an array ebins[0..nbins-1] containing the expected
        // numbers of events, and given the number of constraints knstrn
        // (normally one), this routine returns (trivially) the number of
        // degrees of freedom df, and (nontrivially) the chi-square chsq and
        // the p-value prob. A small value of prob indicates a significant
        // difference between the distributions bins and ebins. Note that bins
        // and ebins are both double arrays, although bins will normally
        // contain integer values.
        Gamma gam = new Gamma();
        int j, nbins = bins.length;
        double temp;
        df.$(nbins - knstrn);
        chsq.$(0.0);
        for (j = 0; j < nbins; j++) {
            if (ebins[j] < 0.0 || (ebins[j] == 0. && bins[j] > 0.))
                throw new NRException("Bad expected number in chsone");
            if (ebins[j] == 0.0 && bins[j] == 0.0) {
                df.$(df.$() - 1); // No data means one less degree of free
            } else { // dom.
                temp = bins[j] - ebins[j];
                chsq.$(chsq.$() + temp * temp / ebins[j]);
            }
        }
        prob.$(gam.gammq(0.5 * df.$(), 0.5 * chsq.$())); // Chi-square
                                                                     // probability
                                                                     // function.
                                                                     // See
                                                                     // 6.2.
    }

    public static void chstwo(final double[] bins1, final double[] bins2, 
            final $double df, final $double chsq, final $double prob) throws NRException {
        chstwo(bins1, bins2, df, chsq, prob, 1);
    }

    public static void chstwo(final double[] bins1, final double[] bins2, 
            final $double df, final $double chsq, final $double prob, final int knstrn) throws NRException {
        // Given the arrays bins1[0..nbins-1] and bins2[0..nbins-1], containing
        // two sets of binned data, and given the number of constraints knstrn
        // (normally 1 or 0), this routine returns the number of degrees of
        // freedom df, the chi-square chsq, and the p-value prob. A small value
        // of prob indicates a significant difference between the distributions
        // bins1 and bins2. Note that bins1 and bins2 are both double arrays,
        // although they will normally contain integer values.
        Gamma gam = new Gamma();
        int j, nbins = bins1.length;
        double temp;
        df.$(nbins - knstrn);
        chsq.$(0.0);
        for (j = 0; j < nbins; j++)
            if (bins1[j] == 0.0 && bins2[j] == 0.0)
                df.$(df.$() - 1); // No data means one less degree of free
            else { // dom.
                temp = bins1[j] - bins2[j];
                chsq.$(chsq.$() + (temp * temp / (bins1[j] + bins2[j])));
            }
        prob.$(gam.gammq(0.5 * df.$(), 0.5 * chsq.$())); // Chi-square
                                                                     // probability
                                                                     // function.
                                                                     // See
                                                                     // 6.2.
    }

    public static void ksone(final double[] data, Func_Doub_To_Doub func, 
            final $double d, final $double prob) throws NRException {
        // Given an array data[0..n-1], and given a user-supplied function of a
        // single variable func that is a cumulative distribution function
        // ranging from 0 (for smallest values of its argument) to 1 (for
        // largest values of its argument), this routine returns the K–S
        // statistic d and the p-value prob. Small values of prob show that
        // the cumulative distribution function of data is significantly
        // different from func. The array data is modified by being sorted into
        // ascending order.
        int j, n = data.length;
        double dt, en, ff, fn, fo = 0.0;
        KSdist ks = new KSdist();
        sort(data); // If the data are already sorted into ascending
        // order, then this call can be omitted.
        en = n;
        d.$(0.0);
        for (j = 0; j < n; j++) { // Loop over the sorted data points.
            fn = (j + 1) / en; // Data’s c.d.f. after this step.
            ff = func.eval(data[j]); // Compare to the user-supplied function.
            dt = MAX(abs(fo - ff), abs(fn - ff)); // Maximum distance.
            if (dt > d.$())
                d.$(dt);
            fo = fn;
        }
        en = sqrt(en);
        prob.$(ks.qks((en + 0.12 + 0.11 / en) * d.$())); // Compute p-value.
    }

    public static void kstwo(final double[] data1, final double[] data2, 
            final $double d, final $double prob) throws NRException {
        // Given an array data1[0..n1-1], and an array data2[0..n2-1], this
        // routine returns the K–S statistic d and the p-value prob for the
        // null hypothesis that the data sets are drawn from the same
        // distribution. Small values of prob show that the cumulative
        // distribution function of data1 is significantly different from that
        // of data2. The arrays data1 and data2 are modified by being sorted
        // into ascending order.
        int j1 = 0, j2 = 0, n1 = data1.length, n2 = data2.length;
        double d1, d2, dt, en1, en2, en, fn1 = 0.0, fn2 = 0.0;
        KSdist ks = new KSdist();
        sort(data1);
        sort(data2);
        en1 = n1;
        en2 = n2;
        d.$(0.0);
        while (j1 < n1 && j2 < n2) { // If we are not done...
            if ((d1 = data1[j1]) <= (d2 = data2[j2])) // Next step is in data1.
                do
                    fn1 = ++j1 / en1;
                while (j1 < n1 && d1 == data1[j1]);
            if (d2 <= d1) // Next step is in data2.
                do
                    fn2 = ++j2 / en2;
                while (j2 < n2 && d2 == data2[j2]);
            if ((dt = abs(fn2 - fn1)) > d.$())
                d.$(dt);
        }
        en = sqrt(en1 * en2 / (en1 + en2));
        prob.$(ks.qks((en + 0.12 + 0.11 / en) * d.$())); // Compute p-value.
    }

    public static void cntab(final double[][] nn, 
            final $double chisq, final $double df, final $double prob, 
            final $double cramrv, final $double ccc) throws NRException {
        // Given a two-dimensional contingency table in the form of an array
        // nn[0..ni-1][0..nj-1] of integers, this routine returns the chi-square
        // chisq, the number of degrees of freedom df, the pvalue prob (small
        // values indicating a significant association), and two measures of
        // association, Cramer’s V (cramrv) and the contingency coefficient C
        // (ccc).
        final double TINY = 1.0e-30; // A small number.
        Gamma gam = new Gamma();
        int i, j, nnj, nni, minij, ni = nn.length, nj = nn.length;
        double sum = 0.0, expctd, temp;
        final double[] sumi = doub_vec(ni), sumj = doub_vec(nj);
        nni = ni; // Number of rows...
        nnj = nj; // ...and columns.
        for (i = 0; i < ni; i++) { // Get the row totals.
            sumi[i] = 0.0;
            for (j = 0; j < nj; j++) {
                sumi[i] += nn[i][j];
                sum += nn[i][j];
            }
            if (sumi[i] == 0.0)
                --nni; // Eliminate any zero rows by reducing the num
        } // ber.
        for (j = 0; j < nj; j++) { // Get the column totals.
            sumj[j] = 0.0;
            for (i = 0; i < ni; i++)
                sumj[j] += nn[i][j];
            if (sumj[j] == 0.0)
                --nnj; // Eliminate any zero columns.
        }
        df.$(nni * nnj - nni - nnj + 1); // Corrected number of degrees of
                                               // freedom.
        chisq.$(0.0);
        for (i = 0; i < ni; i++) { // Do the chi-square sum.
            for (j = 0; j < nj; j++) {
                expctd = sumj[j] * sumi[i] / sum;
                temp = nn[i][j] - expctd;
                chisq.$(chisq.$() + temp * temp / (expctd + TINY)); // Here TINY
                                                               // guarantees
                                                               // that any
                // eliminated row or column will not contribute to the sum.
            }
        }
        prob.$(gam.gammq(0.5 * df.$(), 0.5 * chisq.$())); // Chi-square
                                                                      // probability
                                                                      // function.
        minij = nni < nnj ? nni - 1 : nnj - 1;
        cramrv.$(sqrt(chisq.$() / (sum * minij)));
        ccc.$(sqrt(chisq.$() / (chisq.$() + sum)));
    }

    public static void pearsn(final double[] x, final double[] y, 
            final $double r, final $double prob, final $double z) 
                    throws NRException {
        // Given two arrays x[0..n-1] and y[0..n-1], this routine computes
        // their correlation coefficient r (returned as r), the p-value at
        // which the null hypothesis of zero correlation is disproved (prob
        // whose small value indicates a significant correlation), and
        // Fisher’s z (returned as z), whose value can be used in further
        // statistical tests as described above.
        final double TINY = 1.0e-20; // Will regularize the unusual case of
        Beta beta = new Beta(); // complete correlation.
        int j, n = x.length;
        double yt, xt, t, df;
        double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;
        for (j = 0; j < n; j++) { // Find the means.
            ax += x[j];
            ay += y[j];
        }
        ax /= n;
        ay /= n;
        for (j = 0; j < n; j++) { // Compute the correlation coefficient.
            xt = x[j] - ax;
            yt = y[j] - ay;
            sxx += xt * xt;
            syy += yt * yt;
            sxy += xt * yt;
        }
        r.$(sxy / (sqrt(sxx * syy) + TINY));
        z.$(0.5 * log((1.0 + r.$() + TINY) / (1.0 - r.$() + TINY))); // Fisher’s
                                                                     // z transformation.
        df = n - 2;
        t = r.$() * sqrt(df / ((1.0 - r.$() + TINY) * (1.0 + r.$() + TINY))); // Equation
                                                                                       // (14.5.5).
        prob.$(beta.betai(0.5 * df, 0.5, df / (df + t * t))); // Student’s
                                                                    // t
                                                                    // probability.
        // prob=erfcc(abs(z*sqrt(n-1.0))/1.4142136);
        // For large n, this easier computation of prob, using the short routine
        // erfcc, would give approximately the same value.
    }

    public static void spear(final double[] data1, final double[] data2, 
            final $double d, final $double zd, final $double probd, 
            final $double rs, final $double probrs) 
                    throws NRException {
        // Given two data arrays, data1[0..n-1] and data2[0..n-1], this routine
        // returns their sum squared difference of ranks as D, the number of
        // standard deviations by which D deviates from its null-hypothesis
        // expected value as zd, the two-sided p-value of this deviation as
        // probd, Spearman’s rank correlation rs as rs, and the two-sided
        // p-value of its deviation from zero as probrs. The external routines
        // crank (below) and sort2 (8.2) are used. A small value of either
        // probd or probrs indicates a significant correlation (rs positive) or
        // anticorrelation (rs negative).
        Beta bet = new Beta();
        int j, n = data1.length;
        double vard, t, fac, en3n, en, df, aved;
        $double sf = $(0.0);
        $double sg = $(0.0);
        final double[] wksp1 = doub_vec(n), wksp2 = doub_vec(n);
        for (j = 0; j < n; j++) {
            wksp1[j] = data1[j];
            wksp2[j] = data2[j];
        }
        sort2(wksp1, wksp2); // Sort each of the data arrays, and convert the
                             // entries
        // to ranks. The values sf and sg return the sums
        // P
        // .f 3
        // k  fk/ and
        // P
        // .g3m
        // gm/,
        // respectively.
        crank(wksp1, sf);
        sort2(wksp2, wksp1);
        crank(wksp2, sg);
        d.$(0.0);
        for (j = 0; j < n; j++)
            // Sum the squared difference of ranks.
            d.$(d.$() + SQR(wksp1[j] - wksp2[j]));
        en = n;
        en3n = en * en * en - en;
        aved = en3n / 6.0 - (sf.$() + sg.$()) / 12.0; // Expectation value
                                                            // of D,
        fac = (1.0 - sf.$() / en3n) * (1.0 - sg.$() / en3n);
        vard = ((en - 1.0) * en * en * SQR(en + 1.0) / 36.0) * fac; // and
                                                                    // variance
                                                                    // of D give
        zd.$((d.$() - aved) / sqrt(vard)); // number of standard
                                                    // deviaprobd=
        erfcc(abs(zd.$()) / 1.4142136); // tions and p-value.
        rs.$((1.0 - (6.0 / en3n) * (d.$() + (sf.$() + sg.$()) / 12.0)) / sqrt(fac)); // Rank
                                                                                                    // correlation
                                                                                                    // coefficient,
        fac = (rs.$() + 1.0) * (1.0 - rs.$());
        if (fac > 0.0) {
            t = rs.$() * sqrt((en - 2.0) / fac); // and its t-value,
            df = en - 2.0;
            probrs.$(bet.betai(0.5 * df, 0.5, df / (df + t * t))); // give
                                                                         // its
                                                                         // p-value.
        } else
            probrs.$(0.0);
    }

    public static void crank(final double[] w, final $double s) {
        // Given a sorted array w[0..n-1], replaces the elements by their rank,
        // including midranking of ties, and returns as s the sum of f 3  f ,
        // where f is the number of elements in each tie.
        int j = 1, ji, jt, n = w.length;
        double t, rank;
        s.$(0.0);
        while (j < n) {
            if (w[j] != w[j - 1]) { // Not a tie.
                w[j - 1] = j;
                ++j;
            } else { // A tie:
                for (jt = j + 1; jt <= n && w[jt - 1] == w[j - 1]; jt++)
                    ; // How far does it go?
                rank = 0.5 * (j + jt - 1); // This is the mean rank of the tie,
                for (ji = j; ji <= (jt - 1); ji++)
                    // so enter it into all the tied entries,
                    w[ji - 1] = rank;
                t = jt - j;
                s.$(s.$() + (t * t * t - t)); // and update s.
                j = jt;
            }
        }
        if (j == n)
            w[n - 1] = n; // If the last element was not tied, this is its
    } // rank.

    public static void kendl1(final double[] data1, final double[] data2, 
            final $double tau, final $double z, final $double prob) {
        // Given data arrays data1[0..n-1] and data2[0..n-1], this program
        // returns Kendall’s  as tau, its number of standard deviations from
        // zero as z, and its two-sided p-value as prob. Small values of prob
        // indicate a significant correlation (tau positive) or anticorrelation
        // (tau negative).
        int is = 0, j, k, n2 = 0, n1 = 0, n = data1.length;
        double svar, aa, a2, a1;
        for (j = 0; j < n - 1; j++) { // Loop over first member of pair,
            for (k = j + 1; k < n; k++) { // and second member.
                a1 = data1[j] - data1[k];
                a2 = data2[j] - data2[k];
                aa = a1 * a2;
                if (aa != 0.0) { // Neither array has a tie.
                    ++n1;
                    ++n2;
                    if (aa > 0.) {
                        ++is;
                    } else {
                        --is;
                    }
                    ;
                } else { // One or both arrays have ties.
                    if (a1 != 0.0)
                        ++n1; // An “extra x” event.
                    if (a2 != 0.0)
                        ++n2; // An “extra y” event.
                }
            }
        }
        tau.$(is / (sqrt(Doub(n1)) * sqrt(Doub(n2)))); // Equation
                                                             // (14.6.8).
        svar = (4.0 * n + 10.0) / (9.0 * n * (n - 1.0)); // Equation (14.6.9).
        z.$(tau.$() / sqrt(svar));
        prob.$(erfcc(abs(z.$()) / 1.4142136)); // p-value.
    }

    public static void kendl2(final double[][] tab, 
            final $double tau, final $double z, final $double prob) {
        // Given a two-dimensional table tab[0..i-1][0..j-1], such that
        // tab[k][l] contains the number of events falling in bin k of one
        // variable and bin l of another, this program returns Kendall’s  as
        // tau, its number of standard deviations from zero as z, and its
        // two-sided p-value as prob. Small values of prob indicate a
        // significant correlation (tau positive) or anticorrelation (tau
        // negative) between the two variables. Although tab is a double
        // array, it will normally contain integral values.
        int k, l, nn, mm, m2, m1, lj, li, kj, ki, i = nrows(tab), j = ncols(tab);
        double svar, s = 0.0, points, pairs, en2 = 0.0, en1 = 0.0;
        nn = i * j; // Total number of entries in contingency table.
        points = tab[i - 1][j - 1];
        for (k = 0; k <= nn - 2; k++) { // Loop over entries in table,
            ki = (k / j); // decoding a row,
            kj = k - j * ki; // and a column.
            points += tab[ki][kj]; // Increment the total count of events.
            for (l = k + 1; l <= nn - 1; l++) { // Loop over other member of the
                                                // pair,
                li = l / j; // decoding its row
                lj = l - j * li; // and column.
                mm = (m1 = li - ki) * (m2 = lj - kj);
                pairs = tab[ki][kj] * tab[li][lj];
                if (mm != 0) { // Not a tie.
                    en1 += pairs;
                    en2 += pairs;
                    s += (mm > 0 ? pairs : -pairs); // Concordant, or
                                                    // discordant.
                } else {
                    if (m1 != 0)
                        en1 += pairs;
                    if (m2 != 0)
                        en2 += pairs;
                }
            }
        }
        tau.$(s / sqrt(en1 * en2));
        svar = (4.0 * points + 10.0) / (9.0 * points * (points - 1.0));
        z.$(tau.$() / sqrt(svar));
        prob.$(erfcc(abs(z.$()) / 1.4142136));
    }

    public static void ks2d1s(final double[] x1, final double[] y1,
            final Func_Doub_Doub_DoubRef_DoubRef_DoubRef_DoubRef_DoubRef_DoubRef quadvl, 
            final $double d1, final $double prob) throws NRException {
        // Two-dimensional Kolmogorov-Smirnov test of one sample against a
        // model. Given the x and y coordinates of n1 data points in arrays
        // x1[0..n1-1] and y1[0..n1-1], and given a usersupplied function
        // quadvl that exemplifies the model, this routine returns the
        // two-dimensional K-S statistic as d1, and its p-value as prob. Small
        // values of prob show that the sample is significantly different from
        // the model. Note that the test is slightly distribution-dependent,
        // so prob is only an estimate.
        int j, n1 = x1.length;
        double rr, sqen;
        final $double dum = $(0.0), dumm = $(0.0);
        final $double r1 = $(0.0);
        final $double ga = $(0.0), gb = $(0.0), gc = $(0.0), gd = $(0.0);
        final $double fa = $(0.0), fb = $(0.0), fc = $(0.0), fd = $(0.0);
        KSdist ks = new KSdist();
        d1.$(0.0);
        for (j = 0; j < n1; j++) { // Loop over the data points.
            quadct(x1[j], y1[j], x1, y1, fa, fb, fc, fd);
            quadvl(x1[j], y1[j], ga, gb, gc, gd);
            if (fa.$() > ga.$())
                fa.$(fa.$() + 1.0 / n1);
            if (fb.$() > gb.$())
                fb.$(fb.$() + 1.0 / n1);
            if (fc.$() > gc.$())
                fc.$(fc.$() + 1.0 / n1);
            if (fd.$() > gd.$())
                fd.$(fd.$() + 1.0 / n1);
            d1.$(MAX(d1.$(), abs(fa.$() - ga.$())));
            d1.$(MAX(d1.$(), abs(fb.$() - gb.$())));
            d1.$(MAX(d1.$(), abs(fc.$() - gc.$())));
            d1.$(MAX(d1.$(), abs(fd.$() - gd.$())));
            // For both the sample and the model, the distribution is integrated
            // in each of four quadrants, and the maximum difference is saved.
        }
        pearsn(x1, y1, r1, dum, dumm); // Get the linear correlation
                                                   // coefficient r1.
        sqen = sqrt(Doub(n1));
        rr = sqrt(1.0 - r1.$() * r1.$());
        // Estimate the probability using the K-S probability function.
        prob.$(ks.qks(d1.$() * sqen / (1.0 + rr * (0.25 - 0.75 / sqen))));
    }

    public static void quadct(final double x, final double y, final double[] xx, final double[] yy,
            final $double fa, final $double fb, final $double fc, final $double fd) {
        // Given an origin .x; y/, and an array of nn points with coordinates
        // xx[0..nn-1] and yy[0..nn-1], count how many of them are in each
        // quadrant around the origin, and return the normalized fractions.
        // Quadrants are labeled alphabetically, counterclockwise from the
        // upper right. Used by ks2d1s and ks2d2s.
        int k, na, nb, nc, nd, nn = xx.length;
        double ff;
        na = nb = nc = nd = 0;
        for (k = 0; k < nn; k++) {
            if (yy[k] == y && xx[k] == x)
                continue;
            if (yy[k] > y) {
                if (xx[k] > x) {
                    ++na;
                } else {
                    ++nb;
                }
                ;
            } else {
                if (xx[k] > x) {
                    ++nd;
                } else {
                    ++nc;
                }
                ;
            }
        }
        ff = 1.0 / nn;
        fa.$(ff * na);
        fb.$(ff * nb);
        fc.$(ff * nc);
        fd.$(ff * nd);
    }

    public static void quadvl(final double x, final double y, 
            final $double fa, final $double fb,
            final $double fc, final $double fd) {
        // This is a sample of a user-supplied routine to be used with ks2d1s.
        // In this case, the model distribution is uniform inside the
        // square 1 < x < 1, 1 < y < 1. In general, this routine should
        // return,
        // for any point .x; y/, the fraction of the total distribution in
        // each of the four quadrants around that point. The fractions, fa, fb,
        // fc, and fd, must add up to 1. Quadrants are alphabetical,
        // counterclockwise from the upper right.
        double qa, qb, qc, qd;
        qa = MIN(2.0, MAX(0.0, 1.0 - x));
        qb = MIN(2.0, MAX(0.0, 1.0 - y));
        qc = MIN(2.0, MAX(0.0, x + 1.0));
        qd = MIN(2.0, MAX(0.0, y + 1.0));
        fa.$(0.25 * qa * qb);
        fb.$(0.25 * qb * qc);
        fc.$(0.25 * qc * qd);
        fd.$(0.25 * qd * qa);
    }

    // The routine ks2d2s is the two-sample case of the two-dimensional K–S
    // test. It also calls quadct, pearsn, and KSdist::qks. Being a two-sample
    // test, it does not need an analytic model.

    public static void ks2d2s(final double[] x1, final double[] y1, final double[] x2, final double[] y2,
            final $double d, final $double prob) throws NRException {
        // Two-dimensional Kolmogorov-Smirnov test on two samples. Given the x
        // and y coordinates of the first sample as n1 values in arrays
        // x1[0..n1-1] and y1[0..n1-1], and likewise for the second sample,
        // n2 values in arrays x2 and y2, this routine returns the
        // two-dimensional, twosample K-S statistic as d, and its p-value as
        // prob. Small values of prob show that the two samples are
        // significantly different. Note that the test is slightly
        // distribution-dependent, so prob is only an estimate.
        int j, n1 = x1.length, n2 = x2.length;
        double d1, d2, rr, sqen;
        final $double dum = $(0.0), dumm = $(0.0), r1 = $(0.0), r2 = $(0.0);
        final $double fa = $(0.0), fb = $(0.0), fc = $(0.0), fd = $(0.0), ga = $(0.0);
        final $double gb = $(0.0), gc = $(0.0), gd = $(0.0);
        KSdist ks = new KSdist();
        d1 = 0.0;
        for (j = 0; j < n1; j++) { // First, use points in the first sample as
                                   // origins.
            quadct(x1[j], y1[j], x1, y1, fa, fb, fc, fd);
            quadct(x1[j], y1[j], x2, y2, ga, gb, gc, gd);
            if (fa.$() > ga.$())
                fa.$(fa.$() + 1.0 / n1);
            if (fb.$() > gb.$())
                fb.$(fb.$() + 1.0 / n1);
            if (fc.$() > gc.$())
                fc.$(fc.$() + 1.0 / n1);
            if (fd.$() > gd.$())
                fd.$(fd.$() + 1.0 / n1);
            d1 = MAX(d1, abs(fa.$() - ga.$()));
            d1 = MAX(d1, abs(fb.$() - gb.$()));
            d1 = MAX(d1, abs(fc.$() - gc.$()));
            d1 = MAX(d1, abs(fd.$() - gd.$()));
        }
        d2 = 0.0;
        for (j = 0; j < n2; j++) { // Then, use points in the second sample as
            quadct(x2[j], y2[j], x1, y1, fa, fb, fc, fd); // origins.
            quadct(x2[j], y2[j], x2, y2, ga, gb, gc, gd);
            if (ga.$() > fa.$())
                ga.$(ga.$() + 1.0 / n1);
            if (gb.$() > fb.$())
                gb.$(gb.$() + 1.0 / n1);
            if (gc.$() > fc.$())
                gc.$(gc.$() + 1.0 / n1);
            if (gd.$() > fd.$())
                gd.$(gd.$() + 1.0 / n1);
            d2 = MAX(d2, abs(fa.$() - ga.$()));
            d2 = MAX(d2, abs(fb.$() - gb.$()));
            d2 = MAX(d2, abs(fc.$() - gc.$()));
            d2 = MAX(d2, abs(fd.$() - gd.$()));
        }
        d.$(0.5 * (d1 + d2)); // Average the K-S statistics.
        sqen = sqrt(n1 * n2 / Doub(n1 + n2));
        pearsn(x1, y1, r1, dum, dumm); // Get the linear correlation
                                                   // coefficient for each
        pearsn(x2, y2, r2, dum, dumm); // sample.
        rr = sqrt(1.0 - 0.5 * (r1.$() * r1.$() + r2.$() * r2.$()));
        // Estimate the probability using the K-S probability function.
        prob.$(ks.qks(d.$() * sqen / (1.0 + rr * (0.25 - 0.75 / sqen))));
    }

    public static void savgol(final double[] c, final int np, final int nl, final int nr, final int ld, final int m)
            throws NRException {
        // Returns in c[0..np-1], in wraparound order (N.B.!) consistent with
        // the argument respns in routine convlv, a set of Savitzky-Golay
        // filter coefficients. nl is the number of leftward (past) data points
        // used, while nr is the number of rightward (future) data points,
        // making the total number of data points used nlCnrC1. ld is the order
        // of the derivative desired (e.g., ld D 0 for smoothed function. For
        // the derivative of order k, you must multiply the array c by kŠ.) m
        // is the order of the smoothing polynomial, also equal to the highest
        // conserved moment; usual values are m D 2 or m D 4.
        int j, k, imj, ipj, kk, mm;
        double fac, sum;
        if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
            throw new NRException("bad args in savgol");
        @SuppressWarnings("unused")
        final int[] indx = int_vec(m + 1);
        final double[][] a = doub_mat(m + 1, m + 1);
        final double[] b = doub_vec(m + 1);
        for (ipj = 0; ipj <= (m << 1); ipj++) { // Set up the normal equations
                                                // of the desired
            sum = ((ipj != 0) ? 0.0 : 1.0); // least-squares fit.
            for (k = 1; k <= nr; k++)
                sum += pow(Doub(k), Doub(ipj));
            for (k = 1; k <= nl; k++)
                sum += pow(Doub(-k), Doub(ipj));
            mm = MIN(ipj, 2 * m - ipj);
            for (imj = -mm; imj <= mm; imj += 2)
                a[(ipj + imj) / 2][(ipj - imj) / 2] = sum;
        }
        LUdcmp alud = new LUdcmp(a); // Solve them: LU decomposition.
        for (j = 0; j < m + 1; j++)
            b[j] = 0.0;
        b[ld] = 1.0;
        // Right-hand side vector is unit vector, depending on which derivative
        // we want.
        alud.solve(b, b); // Get one row of the inverse matrix.
        for (kk = 0; kk < np; kk++)
            c[kk] = 0.0; // Zero the output array (it may be bigger than
        for (k = -nl; k <= nr; k++) { // number of coefficients).
            sum = b[0]; // Each Savitzky-Golay coefficient is the dot
            // product of powers of an integer with the inverse matrix row.
            fac = 1.0;
            for (mm = 1; mm <= m; mm++)
                sum += b[mm] * (fac *= k);
            kk = (np - k) % np; // Store in wraparound order.
            c[kk] = sum;
        }
    }

}
