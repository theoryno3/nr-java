package com.snuggy.nr.chapter10;

import static java.lang.Math.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Amoeba {

	// Multidimensional minimization by the downhill simplex method of Nelder
	// and Mead.
	private final double ftol;
	private int nfunc; // The number of function evaluations.
	private int mpts;
	private int ndim;
	@SuppressWarnings("unused")
	private double fmin; // Function value at the minimum.
	private double[] y; // Function values at the vertices of the simplex.
	private double[][] p; // Current simplex.

	public Amoeba(final double ftoll) {
		ftol = (ftoll);
	}

	// The Constructor argument ftoll is the fractional convergence tolerance
	// to be achieved in the function value (n.b.!).

	public <T extends Func_DoubVec_To_Doub> double[] minimize(
			final double[] point, final double del, final T func)
			throws NRException {
		// Multidimensional minimization of the function or functor func(x),
		// where x[0..ndim-1] is a vector in ndim dimensions, by the downhill
		// simplex method of Nelder and Mead. The initial simplex is specified
		// as in equation (10.5.1) by a point[0..ndim-1] and a constant
		// displacement del along each coordinate direction. Returned is the
		// location of the minimum.
		double[] dels = doub_vec(point.length, del);
		return minimize(point, dels, func);
	}

	public <T extends Func_DoubVec_To_Doub> double[] minimize(
			final double[] point, final double[] dels, final T func)
			throws NRException {
		// Alternative interface that takes different displacements
		// dels[0..ndim-1] in different directions for the initial simplex.
		int ndim = point.length;
		double[][] pp = doub_mat(ndim + 1, ndim);
		for (int i = 0; i < ndim + 1; i++) {
			for (int j = 0; j < ndim; j++)
				pp[i][j] = point[j];
			if (i != 0)
				pp[i][i - 1] += dels[i - 1];
		}
		return minimize(pp, func);
	}

	public <T extends Func_DoubVec_To_Doub> double[] minimize(
			final double[][] pp, final T func) throws NRException {
		// Most general interface: initial simplex specified by the matrix
		// pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors
		// that are the vertices of the starting simplex.
		final int NMAX = 5000; // Maximum allowed number of function evalua
		final double TINY = 1.0e-10; // tions.
		int ihi, ilo, inhi;
		mpts = nrows(pp);
		ndim = ncols(pp);
		double[] psum = doub_vec(ndim), pmin = doub_vec(ndim), x = doub_vec(ndim);
		p = pp;
		// y.resize(mpts);
		y = new double[mpts];
		for (int i = 0; i < mpts; i++) {
			for (int j = 0; j < ndim; j++)
				x[j] = p[i][j];
			y[i] = func.eval(x);
		}
		nfunc = 0;
		get_psum(p, psum);
		for (;;) {
			ilo = 0;
			// First we must determine which point is the highest (worst),
			// next-highest, and lowest (best), by looping over the points
			// in the simplex.
			// ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			if (y[0] > y[1]) {
				inhi = 1;
				ihi = 0;
			} else {
				inhi = 0;
				ihi = 1;
			}
			for (int i = 0; i < mpts; i++) {
				if (y[i] <= y[ilo])
					ilo = i;
				if (y[i] > y[ihi]) {
					inhi = ihi;
					ihi = i;
				} else if (y[i] > y[inhi] && i != ihi)
					inhi = i;
			}
			double rtol = 2.0 * abs(y[ihi] - y[ilo])
					/ (abs(y[ihi]) + abs(y[ilo]) + TINY);
			// Compute the fractional range from highest to lowest and return if
			// satisfactory.
			// I think p[ilo][*] is the current argmin here.
			for (int i = 0; i < ndim; i++)
				System.out.printf("%f ", p[ilo][i]);
			System.out.printf("\n");
			if (rtol < ftol) { // If returning, put best point and value in slot
								// 0.
				SWAP(y, 0, ilo);
				for (int i = 0; i < ndim; i++) {
					SWAP(p, 0, i, ilo, i);
					pmin[i] = p[0][i];
				}
				fmin = y[0];
				return pmin;
			}
			if (nfunc >= NMAX)
				throw new NRException("NMAX exceeded");
			nfunc += 2;
			// Begin a new iteration. First extrapolate by a factor 1 through
			// the
			// face of the simplex across from the high point, i.e., reflect the
			// simplex from the high point.
			double ytry = amotry(p, y, psum, ihi, -1.0, func);
			if (ytry <= y[ilo])
				// Gives a result better than the best point, so try an
				// additional
				// extrapolation by a factor 2.
				ytry = amotry(p, y, psum, ihi, 2.0, func);
			else if (ytry >= y[inhi]) {
				// The reflected point is worse than the second-highest, so look
				// for
				// an intermediate lower point, i.e., do a one-dimensional
				// contraction.
				double ysave = y[ihi];
				ytry = amotry(p, y, psum, ihi, 0.5, func);
				if (ytry >= ysave) { // Can't seem to get rid of that high
										// point.
					// Better contract around the lowest (best) point.
					for (int i = 0; i < mpts; i++) {
						if (i != ilo) {
							for (int j = 0; j < ndim; j++)
								p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							y[i] = func.eval(psum);
						}
					}
					nfunc += ndim; // Keep track of function evaluations.
					get_psum(p, psum); // Recompute psum.
				}
			} else
				--nfunc; // Correct the evaluation count.
		} // Go back for the test of doneness and the next
	} // iteration.

	public void get_psum(final double[][] p, final double[] psum) {
		// Utility function.
		for (int j = 0; j < ndim; j++) {
			double sum = 0.0;
			for (int i = 0; i < mpts; i++)
				sum += p[i][j];
			psum[j] = sum;
		}
	}

	public <T extends Func_DoubVec_To_Doub> double amotry(final double[][] p,
			final double[] y, final double[] psum, final int ihi,
			final double fac, final T func) throws NRException {
		// Helper function: Extrapolates by a factor fac through the face of the
		// simplex across from the high point, tries it, and replaces the high
		// point if the new point is better.
		double[] ptry = doub_vec(ndim);
		double fac1 = (1.0 - fac) / ndim;
		double fac2 = fac1 - fac;
		for (int j = 0; j < ndim; j++)
			ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
		double ytry = func.eval(ptry); // Evaluate the function at the trial
										// point.
		if (ytry < y[ihi]) { // If it's better than the highest, then replace
								// the
			y[ihi] = ytry; // highest.
			for (int j = 0; j < ndim; j++) {
				psum[j] += ptry[j] - p[ihi][j];
				p[ihi][j] = ptry[j];
			}
		}
		return ytry;
	}
}
