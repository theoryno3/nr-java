package com.snuggy.nr.chapter02;

import static com.snuggy.nr.util.Static.*;

import java.util.*;

public class ADAT {

    final private NRsparseMat a, at; // Store references to A and AT .
    private NRsparseMat adat; // This will hold ADAT .

    // ADAT(const NRsparseMat &A,const NRsparseMat &AT);
    // Allocates compressed column storage for AAT , where A and AT are input
    // in compressed column format, and fills in values of col_ptr and row_ind.
    // Each column must be in sorted order in input matrices. Matrix is output
    // with each column sorted.

    // void updateD(const VecDoub &D); Computes ADAT,
    // where D is a diagonal matrix. This function can be called repeatedly to
    // update ADAT for fixed A.

    // NRsparseMat &ref(); Returns reference to adat,
    // which holds ADAT .

    // ~ADAT();

    public ADAT(final NRsparseMat A, final NRsparseMat AT) {
        a = new NRsparseMat(A);
        at = new NRsparseMat(AT);
        int h, i, j, k, l, nvals, m = AT.ncols();
        int[] done = int_arr(m);
        for (i = 0; i < m; i++)
            // Initialize to not done.
            done[i] = -1;
        nvals = 0; // First pass determines number of nonzeros.
        for (j = 0; j < m; j++) { // Outer loop over columns of AT in AAT .
            for (i = AT.col_ptr()[j]; i < AT.col_ptr()[j + 1]; i++) {
                k = AT.row_ind()[i]; // AT kj ¤ 0. Find column k in first
                                     // matrix, A.
                for (l = A.col_ptr()[k]; l < A.col_ptr()[k + 1]; l++) {
                    h = A.row_ind()[l]; // Ahl ¤ 0.
                    if (done[h] != j) { // Test if contribution already
                                        // included.
                        done[h] = j;
                        nvals++;
                    }
                }
            }
        }
        adat = new NRsparseMat(m, m, nvals); // Allocate storage for ADAT.
        for (i = 0; i < m; i++)
            // Re-initialize.
            done[i] = -1;
        nvals = 0;
        // Second pass: Determine columns of adat. Code is identical to first
        // pass except adat->col_ptr and adat->row_ind get assigned at
        // appropriate places.
        for (j = 0; j < m; j++) {
            adat.col_ptr()[j] = nvals;
            for (i = AT.col_ptr()[j]; i < AT.col_ptr()[j + 1]; i++) {
                k = AT.row_ind()[i];
                for (l = A.col_ptr()[k]; l < A.col_ptr()[k + 1]; l++) {
                    h = A.row_ind()[l];
                    if (done[h] != j) {
                        done[h] = j;
                        adat.row_ind()[nvals] = h;
                        nvals++;
                    }
                }
            }
        }
        adat.col_ptr()[m] = nvals; // Set last value.
        for (j = 0; j < m; j++) { // Sort columns
            i = adat.col_ptr()[j];
            int size = adat.col_ptr()[j + 1] - i;
            if (size > 1) {
                // VecInt col(size,&adat->row_ind[i]);
                int[] col = int_arr(size, adat.row_ind(), i);
                Arrays.sort(col);
                for (k = 0; k < size; k++)
                    adat.row_ind()[i + k] = col[k];
            }
        }
    }

    // The next routine, updateD, actually fills in the values in the val array.
    // It can be called repeatedly to update ADAT for fixed A.

    public void updateD(final double[] D) {
        int h, i, j, k, l, m = a.nrows(), n = a.ncols();
        double[] temp = doub_arr(n), temp2 = doub_arr(m);
        for (i = 0; i < m; i++) { // Outer loop over columns of AT .
            for (j = at.col_ptr()[i]; j < at.col_ptr()[i + 1]; j++) {
                k = at.row_ind()[j]; // Scale elements of each column with D and
                temp[k] = at.val()[j] * D[k]; // store in temp.
            }
            for (j = at.col_ptr()[i]; j < at.col_ptr()[i + 1]; j++) {
                // Go down column again.
                k = at.row_ind()[j];
                for (l = a.col_ptr()[k]; l < a.col_ptr()[k + 1]; l++) {
                    // Go down column k in A.
                    h = a.row_ind()[l];
                    temp2[h] += temp[k] * a.val()[l]; // All terms from temp[k]
                                                      // used here.
                }
            }
            for (j = adat.col_ptr()[i]; j < adat.col_ptr()[i + 1]; j++) {
                // Store temp2 in column of answer.
                k = adat.row_ind()[j];
                adat.val()[j] = temp2[k];
                temp2[k] = 0.0; // Restore temp2.
            }
        }
    }

    // The final two functions are simple. The ref routine returns a reference
    // to the matrix ADAT stored in the structure for other routines that may
    // need to work with it. And the destructor releases the storage.

    public NRsparseMat ref() {
        return adat;
    }

    // ADAT::~ADAT() {
    // delete adat;
    // }

}
