package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Static {

    private static final long c1[] = { 0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L };
    private static final long c2[] = { 0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L };

    public static void psdes(final int lword_ref[], final int rword_ref[]) {
        // Pseudo-DES hashing of the 64-bit word (lword,rword). Both 32-bit
        // arguments are returned hashed on all bits.
        final int NITER = 2;
        long i, ia, ib, iswap, itmph = 0, itmpl = 0;
        for (i = 0; i < NITER; i++) {
            // Perform niter iterations of DES logic, using a simpler
            // (noncryptographic) nonlinear function instead of DES’s.
            ia = (iswap = rword_ref[0]) ^ c1[(int) i]; // The bit-rich constants
                                                       // c1 and (below)
            // c2 guarantee lots of nonlinear mixing.
            itmpl = ia & 0xffff;
            itmph = ia >>> 16;
            ib = itmpl * itmpl + ~(itmph * itmph);
            rword_ref[0] = (int) (lword_ref[0] ^ (((ia = (ib >>> 16) | ((ib & 0xffff) << 16)) ^ c2[(int) i]) + itmpl
                    * itmph));
            lword_ref[0] = (int) iswap;
        }
    }

    public static void hashall(final int[] arr) throws NRException {
        // Replace the array arr by a same-sized hash, all of whose bits depend
        // on all of the bits in arr. Uses psdes for the mutual hash of two
        // 32-bit words.
        int m = arr.length, n = m - 1;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;
        // Incredibly, n is now the next power of 2 m.
        int nb = n, nb2 = n >> 1, j, jb;
        if (n < 2)
            throw new NRException("size must be > 1");
        while (nb > 1) {
            for (jb = 0; jb < n - nb + 1; jb += nb)
                for (j = 0; j < nb2; j++)
                    if (jb + j + nb2 < m) {
                        // psdes(arr[jb + j], arr[jb + j + nb2]);
                        int foo_ref[] = int_ref();
                        int bar_ref[] = int_ref();
                        foo_ref[0] = arr[jb + j];
                        bar_ref[0] = arr[jb + j + nb2];
                        psdes(foo_ref, bar_ref);
                        arr[jb + j] = foo_ref[0];
                        arr[jb + j + nb2] = bar_ref[0];
                    }
            nb = nb2;
            nb2 >>= 1;
        }
        nb2 = n >> 1;
        if (m != n)
            for (j = nb2; j < m; j++) {
                // psdes(arr[j], arr[j - nb2]);
                int foo_ref[] = int_ref();
                int bar_ref[] = int_ref();
                foo_ref[0] = arr[j];
                bar_ref[0] = arr[j - nb2];
                psdes(foo_ref, bar_ref);
                arr[j] = foo_ref[0];
                arr[j - nb2] = bar_ref[0];
            }
        // Final mix needed only if m is not a power of 2.
    }

    public static final double[] torusfuncs(final double[] x) {
        // Return the integrands in equation (7.7.5), with  D 1.
        double den = 1.;
        final double[] f = doub_arr(4);
        f[0] = den;
        for (int i = 1; i < 4; i++)
            f[i] = x[i - 1] * den;
        return f;
    }

    public static boolean torusregion(final double[] x) {
        // Return the inequality (7.7.3).
        return SQR(x[2]) + SQR(sqrt(SQR(x[0]) + SQR(x[1])) - 3.) <= 1.;
    }

    public static final double[] torusmap(final double[] s) {
        // Return the mapping from s to z defined by the last equation in
        // (7.7.7), mapping the other coordinates by the identity map.
        final double[] xx = doub_arr(s);
        xx[2] = 0.2 * log(5. * s[2]);
        return xx;
    }

    public static double torusfunc(final double[] x, final double wgt) {
        double den = exp(5. * x[2]);
        if (SQR(x[2]) + SQR(sqrt(SQR(x[0]) + SQR(x[1])) - 3.) <= 1.)
            return den;
        else
            return 0.;
    }

    private static final int MAXBIT = 30, MAXDIM = 6;
    private static int mdeg[] = { 1, 2, 3, 3, 4, 4 };
    private static int in;
    private static final int[] ix = int_arr(MAXDIM);
    // private static NRvector<Uint*> iu(MAXBIT);

    private static IntPointer[] iu;
    private static final int[] ip = { 0, 1, 1, 2, 1, 4 };
    private static final int[] iv = new int[/*MAXDIM*MAXBIT*/] {
         1,  1,  1,  1,  1,  1, 
         3,  1,  3,  3,  1,  1, 
         5,  7,  7,  3,  3,  5, 
         15, 11, 5,  15, 13, 9, 
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,
    };
    private static double fac;
    
    static {
        iu = new IntPointer[MAXBIT];
    }

    public static void sobseq(final int n, final double[] x) throws NRException {
        // When n is negative, internally initializes a set of MAXBIT direction
        // numbers for each of MAXDIM different Sobol’ sequences. When n is
        // positive (but MAXDIM), returns as the vector x[0..n-1] the next
        // values from n of these sequences. (n must not be changed between
        // initializations.)
        int j, k, l;
        int i, im, ipp;

        if (n < 0) { // Initialize, don’t return a vector.
            for (k = 0; k < MAXDIM; k++)
                ix[k] = 0;
            in = 0;
            if (iv[0] != 1)
                return;
            fac = 1.0 / (1 << MAXBIT);
            // for (j=0,k=0;j<MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
            for (j = 0, k = 0; j < MAXBIT; j++, k += MAXDIM)
                iu[j] = new IntPointer(iv, k);
            // To allow both 1D and 2D addressing.
            for (k = 0; k < MAXDIM; k++) {
                // for (j=0;j<mdeg[k];j++) iu[j][k] <<= (MAXBIT-1-j);
                for (j = 0; j < mdeg[k]; j++)
                    iu[j].set(k, iu[j].ref(k) << (MAXBIT - 1 - j));
                // Stored values only require normalization.
                for (j = mdeg[k]; j < MAXBIT; j++) { // Use the recurrence to
                                                     // get other val
                    ipp = ip[k]; // ues.
                    // i=iu[j-mdeg[k]][k];
                    i = iu[j - mdeg[k]].ref(k);
                    i ^= (i >> mdeg[k]);
                    for (l = mdeg[k] - 1; l >= 1; l--) {
                        // if ((ipp & 1) != 0) i ^= iu[j-l][k];
                        if ((ipp & 1) != 0)
                            i ^= iu[j - l].ref(k);
                        ipp >>= 1;
                    }
                    // iu[j][k]=i;
                    iu[j].set(k, i);
                }
            }
        } else { // Calculate the next vector in the se
            im = in++; // quence.
            for (j = 0; j < MAXBIT; j++) { // Find the rightmost zero bit.
                if (!((im & 1) != 0))
                    break;
                im >>= 1;
            }
            if (j >= MAXBIT)
                throw new NRException("MAXBIT too small in sobseq");
            im = j * MAXDIM;
            for (k = 0; k < MIN(n, MAXDIM); k++) { // XOR the appropriate
                                                   // direction number
                // into each component of the vector and convert to a floating
                // number.
                ix[k] ^= iv[im + k];
                x[k] = ix[k] * fac;
            }
        }
    }

    // Best make everything static, allowing restarts.
    private static final int NDMX = 50, MXDIM = 10, vegas_RANSEED = 5330;
    private static final double ALPH = 1.5, TINY = 1.0e-30;
    private static int i, it, j, k, mds, nd, ndo, ng, npg;
    private static double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti;
    private static double tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt;
    private static final int[] ia = int_arr(MXDIM), kg = int_arr(MXDIM);
    private static final double[] dt = doub_arr(MXDIM), dx = doub_arr(MXDIM), r = doub_arr(NDMX), x = doub_arr(MXDIM),
            xin = doub_arr(NDMX);
    private static final double[][] d = doub_mat(NDMX, MXDIM), di = doub_mat(NDMX, MXDIM), xi = doub_mat(MXDIM, NDMX);
    private static Ran ran_vegas = new Ran(vegas_RANSEED); // Initialize a
                                                           // captive,

    // static random number
    // generator.

    public static void vegas(final double[] regn, final Func_DoubArr_Doub_To_Doub fxn, final int init, final int ncall, final int itmx,
            final int nprn, final double tgral_ref[], final double sd_ref[], final double chi2a_ref[]) {
        // Performs Monte Carlo integration of a user-supplied ndim-dimensional
        // function fxn over a rectangular volume speci ed by regn[0..2*ndim-1],
        // a vector consisting of ndim \lower left" coordinates of the region
        // followed by ndim upper right" coordinates. The integration consists
        // of itmx iterations, each with approximately ncall calls to the
        // function. After each iteration the grid is re ned; more than 5 or 10
        // iterations are rarely useful. The input ag init signals whether this
        // call is a new start or a subsequent call for additional iterations
        // (see comments in the code). The input ag nprn (normally 0) controls
        // the amount of diagnostic output. Returned answers are tgral (the best
        // estimate of the integral), sd (its standard deviation), and chi2a (2
        // per degree of freedom, an indicator of whether consistent results are
        // being obtained). See text for further details.

        int ndim = regn.length / 2;
        if (init <= 0) { // Normal entry. Enter here on a cold start.
            mds = ndo = 1; // Change to mds=0 to disable strati ed sampling,
            for (j = 0; j < ndim; j++)
                xi[j][0] = 1.0; // i.e., use importance sampling only.
        }
        if (init <= 1)
            si = swgt = schi = 0.0;
        // Enter here to inherit the grid from a previous call, but not its
        // answers.
        if (init <= 2) { // Enter here to inherit the previous grid and its
            nd = NDMX; // answers.
            ng = 1;
            if (mds != 0) { // Set up for strati cation.
                ng = Int(pow(ncall / 2.0 + 0.25, 1.0 / ndim));
                mds = 1;
                if ((2 * ng - NDMX) >= 0) {
                    mds = -1;
                    npg = ng / NDMX + 1;
                    nd = ng / npg;
                    ng = npg * nd;
                }
            }
            for (k = 1, i = 0; i < ndim; i++)
                k *= ng;
            npg = MAX(Int(ncall / k), 2);
            calls = Doub(npg) * Doub(k);
            dxg = 1.0 / ng;
            for (dv2g = 1, i = 0; i < ndim; i++)
                dv2g *= dxg;
            dv2g = SQR(calls * dv2g) / npg / npg / (npg - 1.0);
            xnd = nd;
            dxg *= xnd;
            xjac = 1.0 / calls;
            for (j = 0; j < ndim; j++) {
                dx[j] = regn[j + ndim] - regn[j];
                xjac *= dx[j];
            }
            if (nd != ndo) { // Do binning if necessary.
                for (i = 0; i < MAX(nd, ndo); i++)
                    r[i] = 1.0;
                for (j = 0; j < ndim; j++)
                    rebin(ndo / xnd, nd, r, xin, xi, j);
                ndo = nd;
            }
            if (nprn >= 0) {
                // cout << " Input parameters for vegas";
                // cout << " ndim= " << setw(4) << ndim;
                // cout << " ncall= " << setw(8) << calls << endl;
                // cout << setw(34) << " it=" << setw(5) << it;
                // cout << " itmx=" << setw(5) << itmx << endl;
                // cout << setw(34) << " nprn=" << setw(5) << nprn;
                // cout << " ALPH=" << setw(9) << ALPH << endl;
                // cout << setw(34) << " mds=" << setw(5) << mds;
                // cout << " nd=" << setw(5) << nd << endl;
                for (j = 0; j < ndim; j++) {
                    // cout << setw(30) << " x1[" << setw(2) << j;
                    // cout << "]= " << setw(11) << regn[j] << " xu[";
                    // cout << setw(2) << j << "]= ";
                    // cout << setw(11) << regn[j+ndim] << endl;
                }
            }
        }
        for (it = 0; it < itmx; it++) {
            // Main iteration loop. Can enter here (init  3) to do an
            // additional
            // itmx iterations with all other parameters unchanged.
            ti = tsi = 0.0;
            for (j = 0; j < ndim; j++) {
                kg[j] = 1;
                for (i = 0; i < nd; i++)
                    d[i][j] = di[i][j] = 0.0;
            }
            for (;;) {
                fb = f2b = 0.0;
                for (k = 0; k < npg; k++) {
                    wgt = xjac;
                    for (j = 0; j < ndim; j++) {
                        xn = (kg[j] - ran_vegas.doub()) * dxg + 1.0;
                        ia[j] = MAX(MIN(Int(xn), NDMX), 1);
                        if (ia[j] > 1) {
                            xo = xi[j][ia[j] - 1] - xi[j][ia[j] - 2];
                            rc = xi[j][ia[j] - 2] + (xn - ia[j]) * xo;
                        } else {
                            xo = xi[j][ia[j] - 1];
                            rc = (xn - ia[j]) * xo;
                        }
                        x[j] = regn[j] + rc * dx[j];
                        wgt *= xo * xnd;
                    }
                    f = wgt * fxn.eval(x, wgt);
                    f2 = f * f;
                    fb += f;
                    f2b += f2;
                    for (j = 0; j < ndim; j++) {
                        di[ia[j] - 1][j] += f;
                        if (mds >= 0)
                            d[ia[j] - 1][j] += f2;
                    }
                }
                f2b = sqrt(f2b * npg);
                f2b = (f2b - fb) * (f2b + fb);
                if (f2b <= 0.0)
                    f2b = TINY;
                ti += fb;
                tsi += f2b;
                if (mds < 0) { // Use strati ed sampling.
                    for (j = 0; j < ndim; j++)
                        d[ia[j] - 1][j] += f2b;
                }
                for (k = ndim - 1; k >= 0; k--) {
                    kg[k] %= ng;
                    if (++kg[k] != 1)
                        break;
                }
                if (k < 0)
                    break;
            }
            tsi *= dv2g; // Compute nal results for this iteration.
            wgt = 1.0 / tsi;
            si += wgt * ti;
            schi += wgt * ti * ti;
            swgt += wgt;
            tgral_ref[0] = si / swgt;
            chi2a_ref[0] = (schi - si * tgral_ref[0]) / (it + 0.0001);
            if (chi2a_ref[0] < 0.0)
                chi2a_ref[0] = 0.0;
            sd_ref[0] = sqrt(1.0 / swgt);
            tsi = sqrt(tsi);
            if (nprn >= 0) {
                // cout << " iteration no. " << setw(3) << (it+1);
                // cout << " : integral = " << setw(14) << ti;
                // cout << " +/- " << setw(9) << tsi << endl;
                // cout << " all iterations: " << " integral =";
                // cout << setw(14) << tgral << "+-" << setw(9) << sd;
                // cout << " chi**2/IT n =" << setw(9) << chi2a << endl;
                if (nprn != 0) {
                    for (j = 0; j < ndim; j++) {
                        // cout << " DATA FOR axis " << setw(2) << j << endl;
                        // cout << " X delta i X delta i";
                        // cout << " X deltai" << endl;
                        for (i = nprn / 2; i < nd - 2; i += nprn + 2) {
                            // cout << setw(8) << xi[j][i] << setw(12) <<
                            // di[i][j];
                            // cout << setw(12) << xi[j][i+1] << setw(12) <<
                            // di[i+1][j];
                            // cout << setw(12) << xi[j][i+2] << setw(12) <<
                            // di[i+2][j];
                            // cout << endl;
                        }
                    }
                }
            }
            for (j = 0; j < ndim; j++) { // Re ne the grid. Consult references
                                         // to understand
                // the subtlety of this procedure. The re nement
                // is damped, to avoid rapid, destabilizing
                // changes, and also compressed in range
                // by the exponent ALPH.
                xo = d[0][j];
                xn = d[1][j];
                d[0][j] = (xo + xn) / 2.0;
                dt[j] = d[0][j];
                for (i = 2; i < nd; i++) {
                    rc = xo + xn;
                    xo = xn;
                    xn = d[i][j];
                    d[i - 1][j] = (rc + xn) / 3.0;
                    dt[j] += d[i - 1][j];
                }
                d[nd - 1][j] = (xo + xn) / 2.0;
                dt[j] += d[nd - 1][j];
            }
            for (j = 0; j < ndim; j++) {
                rc = 0.0;
                for (i = 0; i < nd; i++) {
                    if (d[i][j] < TINY)
                        d[i][j] = TINY;
                    r[i] = pow((1.0 - d[i][j] / dt[j]) / (log(dt[j]) - log(d[i][j])), ALPH);
                    rc += r[i];
                }
                rebin(rc / xnd, nd, r, xin, xi, j);
            }
        }
    }

    public static void rebin(final double rc, final int nd, final double[] r, final double[] xin, final double[][] xi,
            final int j) {
        // Utility routine used by vegas, to rebin a vector of densities
        // contained in row j of xi into new bins de ned by a vector r.
        int i, k = 0;
        double dr = 0.0, xn = 0.0, xo = 0.0;
        for (i = 0; i < nd - 1; i++) {
            while (rc > dr)
                dr += r[(++k) - 1];
            if (k > 1)
                xo = xi[j][k - 2];
            xn = xi[j][k - 1];
            dr -= rc;
            xin[i] = xn - (xn - xo) * dr / r[k - 1];
        }
        for (i = 0; i < nd - 1; i++)
            xi[j][i] = xin[i];
        xi[j][nd - 1] = 1.0;
    }

    private static int iran = 0;

    public static void miser(final Func_DoubArr_To_Doub func, final double[] regn, final int npts, final double dith,
            final double ave_ref[], final double var_ref[]) {
        // Monte Carlo samples a user-supplied ndim-dimensional function func
        // in a rectangular volume speci ed by regn[0..2*ndim-1], a vector
        // consisting of ndim \lower-left" coordinates of the region followed
        // by ndim upper-right" coordinates. The function is sampled a total of
        // npts times, at locations determined by the method of recursive
        // strati ed sampling. The mean value of the function in the region is
        // returned as ave; an estimate of the statistical uncertainty of ave
        // (square of standard deviation) is returned as var. The input
        // parameter dith should normally be set to zero, but can be set to
        // (e.g.) 0.1 if func's active region falls on the boundary of a
        // power-of-two subdivision of region.
        final int MNPT = 15, MNBS = 60;
        final double PFAC = 0.1, TINY = 1.0e-30, BIG = 1.0e30;
        // Here PFAC is the fraction of remaining function evaluations used at
        // each stage to explore the variance of func. At least MNPT function
        // evaluations are performed in any terminal subregion; a subregion is
        // further bisected only if at least MNBS function evaluations are
        // available. We take MNBS D 4  MNPT.
        int j, jb, n, ndim, npre, nptl, nptr;
        double fracl, fval, rgl, rgm, rgr, s, sigl, siglb, sigr, sigrb;
        double avel_ref[] = doub_ref();
        double varl_ref[] = doub_ref();
        double sum, sumb, summ, summ2;
        ndim = regn.length / 2;
        final double[] pt = doub_arr(ndim);
        if (npts < MNBS) { // Too few points to bisect; do straight
            summ = summ2 = 0.0; // Monte Carlo.
            for (n = 0; n < npts; n++) {
                ranpt(pt, regn);
                fval = func.eval(pt);
                summ += fval;
                summ2 += fval * fval;
            }
            ave_ref[0] = summ / npts;
            var_ref[0] = MAX(TINY, (summ2 - summ * summ / npts) / (npts * npts));
        } else { // Do the preliminary (uniform) sampling.
            final double[] rmid = doub_arr(ndim);
            npre = MAX(Int(npts * PFAC), Int(MNPT));
            final double[] fmaxl = doub_arr(ndim), fmaxr = doub_arr(ndim), fminl = doub_arr(ndim), fminr = doub_arr(ndim);
            for (j = 0; j < ndim; j++) { // Initialize the left and right bounds
                                         // for
                iran = (iran * 2661 + 36979) % 175000; // each dimension.
                s = SIGN(dith, Doub(iran - 87500));
                rmid[j] = (0.5 + s) * regn[j] + (0.5 - s) * regn[ndim + j];
                fminl[j] = fminr[j] = BIG;
                fmaxl[j] = fmaxr[j] = (-BIG);
            }
            for (n = 0; n < npre; n++) { // Loop over the points in the sample.
                ranpt(pt, regn);
                fval = func.eval(pt);
                for (j = 0; j < ndim; j++) { // Find the left and right bounds
                                             // for each
                    if (pt[j] <= rmid[j]) { // dimension.
                        fminl[j] = MIN(fminl[j], fval);
                        fmaxl[j] = MAX(fmaxl[j], fval);
                    } else {
                        fminr[j] = MIN(fminr[j], fval);
                        fmaxr[j] = MAX(fmaxr[j], fval);
                    }
                }
            }
            sumb = BIG; // Choose which dimension jb to bisect.
            jb = -1;
            siglb = sigrb = 1.0;
            for (j = 0; j < ndim; j++) {
                if (fmaxl[j] > fminl[j] && fmaxr[j] > fminr[j]) {
                    sigl = MAX(TINY, pow(fmaxl[j] - fminl[j], 2.0 / 3.0));
                    sigr = MAX(TINY, pow(fmaxr[j] - fminr[j], 2.0 / 3.0));
                    sum = sigl + sigr; // Equation (7.9.24), see text.
                    if (sum <= sumb) {
                        sumb = sum;
                        jb = j;
                        siglb = sigl;
                        sigrb = sigr;
                    }
                }
            }
            if (jb == -1)
                jb = (ndim * iran) / 175000; // MNPT may be too small.
            rgl = regn[jb]; // Apportion the remaining points between
            rgm = rmid[jb]; // left and right.
            rgr = regn[ndim + jb];
            fracl = abs((rgm - rgl) / (rgr - rgl));
            nptl = Int(MNPT + (npts - npre - 2 * MNPT) * fracl * siglb / (fracl * siglb + (1.0 - fracl) * sigrb)); // Equation
                                                                                                                   // (7.9.23).
            nptr = npts - npre - nptl;
            final double[] regn_temp = doub_arr(2 * ndim); // Now allocate and
                                                     // integrate the two sub
            for (j = 0; j < ndim; j++) { // regions.
                regn_temp[j] = regn[j];
                regn_temp[ndim + j] = regn[ndim + j];
            }
            regn_temp[ndim + jb] = rmid[jb];
            miser(func, regn_temp, nptl, dith, avel_ref, varl_ref);
            regn_temp[jb] = rmid[jb]; // Dispatch recursive call; will return
                                      // back
            regn_temp[ndim + jb] = regn[ndim + jb]; // here eventually.
            miser(func, regn_temp, nptr, dith, ave_ref, var_ref);
            ave_ref[0] = fracl * avel_ref[0] + (1 - fracl) * ave_ref[0];
            var_ref[0] = fracl * fracl * varl_ref[0] + (1 - fracl) * (1 - fracl) * var_ref[0];
            // Combine left and right regions by equation (7.9.11) (1st line).
        }
    }

    private static final int ranpt_RANSEED = 5331;
    private static Ran ran = new Ran(ranpt_RANSEED);

    public static void ranpt(final double[] pt, final double[] regn) {
        // Returns a uniformly random point pt in an n-dimensional rectangular
        // region. Used by miser.
        int j, n = pt.length;
        for (j = 0; j < n; j++)
            pt[j] = regn[j] + (regn[n + j] - regn[j]) * ran.doub();
    }

}
