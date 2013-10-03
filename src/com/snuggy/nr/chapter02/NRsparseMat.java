
package com.snuggy.nr.chapter02;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;

public class NRsparseMat implements ByValue<NRsparseMat> {

    // Sparse matrix data structure for compressed column storage.

    private int nrows; // Number of rows.
    private int ncols; // Number of columns.
    private int nvals; // Maximum number of nonzeros.
    private int[] col_ptr; // Pointers to start of columns. Length is ncols+1.
    private int[] row_ind; // Row indices of nonzeros.
    private double[] val; // Array of nonzero values.
    
    @Override
    public NRsparseMat copyOut() {
        NRsparseMat r = new NRsparseMat(this);
	    return r;
    }

    @Override
    public void copyIn(NRsparseMat t) {
	    nrows = t.nrows; 
	    ncols = t.ncols; 
	    nvals = t.nvals; 
	    System.arraycopy(t.col_ptr, 0, col_ptr, 0, t.col_ptr.length);
	    System.arraycopy(t.row_ind, 0, row_ind, 0, t.row_ind.length);
	    System.arraycopy(t.val, 0, val, 0, t.val.length);
    }

    // NRsparseMat(); Default constructor.
    // NRsparseMat(Int m,Int n,Int nnvals); Constructor. Initializes vector to
    // zero.
    // VecDoub ax(const VecDoub &x) const; Multiply A by a vector x[0..ncols-1].
    // VecDoub atx(const VecDoub &x) const; Multiply AT by a vector
    // x[0..nrows-1].
    // NRsparseMat transpose() const; Form AT .

    // The code for the constructors is standard:

    public int nvals() {
        return nvals;
    }

    public int[] row_ind() {
        return row_ind;
    }

    public int[] col_ptr() {
        return col_ptr;
    }

    public double[] val() {
        return val;
    }

    public int nrows() {
        return nrows;
    }

    public int ncols() {
        return ncols;
    }

    /*
    public NRsparseMat() {
        nrows = (0);
        ncols = (0);
        nvals = (0);
        col_ptr = null;
        row_ind = null;
        val = null;
    }
    */

    public NRsparseMat(final int m, final int n, final int nnvals) {
        nrows = (m);
        ncols = (n);
        nvals = (nnvals);
        col_ptr = int_arr(n + 1);
        row_ind = int_arr(nnvals);
        val = doub_arr(nnvals);
    }

    public NRsparseMat(final NRsparseMat mat) {
        nrows = mat.nrows;
        ncols = mat.ncols;
        nvals = mat.nvals;
        col_ptr = int_arr(mat.col_ptr.length);
        System.arraycopy(mat.col_ptr, 0, col_ptr, 0, col_ptr.length);
        row_ind = int_arr(mat.row_ind.length);
        System.arraycopy(mat.row_ind, 0, row_ind, 0, row_ind.length);
        val = doub_arr(mat.val.length);
        System.arraycopy(mat.val, 0, val, 0, val.length);
    }

    // The single most important use of a matrix in compressed column storage
    // mode is to multiply a vector to its right. Don’t implement this by
    // traversing the rows of A, which is extremely inefficient in this storage
    // mode. Here’s the right way to do it:

    public double[] ax(final double[] x) {
        double[] y = doub_arr(nrows);
        for (int j = 0; j < ncols; j++) {
            for (int i = col_ptr[j]; i < col_ptr[j + 1]; i++)
                y[row_ind[i]] += val[i] * x[j];
        }
        return y;
    }

    // Some inefficiency occurs because of the indirect addressing. While there
    // are other storage modes that minimize this, they have their own
    // drawbacks.
    // It is also simple to multiply the transpose of a matrix by a vector to
    // its right, since we just traverse the columns directly. (Indirect
    // addressing is still required.) Note that the transpose matrix is not
    // actually constructed.

    public double[] atx(final double[] x) {
        double[] y = doub_arr(ncols);
        for (int i = 0; i < ncols; i++) {
            y[i] = 0.0;
            for (int j = col_ptr[i]; j < col_ptr[i + 1]; j++)
                y[i] += val[j] * x[row_ind[j]];
        }
        return y;
    }

    // Because the choice of compressed column storage treats rows and columns
    // quite differently, it is rather an involved operation to construct the
    // transpose of a matrix, given the matrix itself in compressed column
    // storage mode. When the operation cannot be avoided, it is

    public $$<NRsparseMat> transpose() {
        int i, j, k, index, m = nrows, n = ncols;
        NRsparseMat at = new NRsparseMat(n, m, nvals); // Initialized to zero.
        // First find the column lengths for AT , i.e. the row lengths of A.
        int[] count = int_arr(m); // Temporary counters for each row of A.
        for (i = 0; i < n; i++)
            for (j = col_ptr[i]; j < col_ptr[i + 1]; j++) {
                k = row_ind[j];
                count[k]++;
            }
        for (j = 0; j < m; j++)
            // Now set at.col_ptr. 0th entry stays 0.
            at.col_ptr[j + 1] = at.col_ptr[j] + count[j];
        for (j = 0; j < m; j++)
            // Reset counters to zero.
            count[j] = 0;
        for (i = 0; i < n; i++)
            // Main loop.
            for (j = col_ptr[i]; j < col_ptr[i + 1]; j++) {
                k = row_ind[j];
                index = at.col_ptr[k] + count[k]; // Element’s position in
                                                  // column of AT .
                at.row_ind[index] = i;
                at.val[index] = val[j];
                count[k]++; // Increment counter for next element in that
            } // column.
        return $$(at);
    }
}
