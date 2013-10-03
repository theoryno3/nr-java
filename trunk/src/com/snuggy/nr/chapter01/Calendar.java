
package com.snuggy.nr.chapter01;

import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Calendar {
    
    public static void flmoon(final int n, final int nph, final int[] jd_ref, final double[] frac_ref) 
            throws NRException {
        // Our routines begin with an introductory comment summarizing their
        // purpose and explaining their calling sequence. This routine
        // calculates the phases of the moon. Given an integer n and a code
        // nph for the phase desired (nph D 0 for new moon, 1 for first
        // quarter, 2 for full, 3 for last quarter), the routine returns the
        // Julian Day Number jd, and the fractional part of a day frac to be
        // added to it, of the nth such phase since January, 1900. Greenwich
        // Mean Time is assumed.

        final double RAD = 3.141592653589793238 / 180.0;
        int i;
        double am, as, c, t, t2, xtra;
        c = n + nph / 4.0; // This is how we comment an individual line.
        t = c / 1236.85;
        t2 = t * t;
        as = 359.2242 + 29.105356 * c; // You aren’t really intended to
                                       // understand
        am = 306.0253 + 385.816918 * c + 0.010730 * t2; // this algorithm, but
                                                        // it does work!
        jd_ref[0] = 2415020 + 28 * n + 7 * nph;
        xtra = 0.75933 + 1.53058868 * c + ((1.178e-4) - (1.55e-7) * t) * t2;
        if (nph == 0 || nph == 2)
            xtra += (0.1734 - 3.93e-4 * t) * sin(RAD * as) - 0.4068 * sin(RAD * am);
        else if (nph == 1 || nph == 3)
            xtra += (0.1721 - 4.0e-4 * t) * sin(RAD * as) - 0.6280 * sin(RAD * am);
        else
            throw new NRException("nph is unknown in flmoon"); // This indicates
                                                             // an error
                                                             // condition.
        i = (int) (xtra >= 0.0 ? floor(xtra) : ceil(xtra - 1.0));
        jd_ref[0] += i;
        frac_ref[0] = xtra - i;
    }

    public static int julday(final int mm, final int id, final int iyyy) throws NRException {
        // In this routine julday returns the Julian Day Number that begins
        // at noon of the calendar date specified by month mm, day id, and
        // year iyyy, all integer variables. Positive year signifies A.D.;
        // negative, B.C. Remember that the year after 1 B.C. was 1 A.D.

        final int IGREG = 15 + 31 * (10 + 12 * 1582); // Gregorian Calendar
                                                      // adopted Oct. 15, 1582.
        int ja, jul, jy = iyyy, jm;
        if (jy == 0)
            throw new NRException("julday: there is no year zero.");
        if (jy < 0)
            ++jy;
        if (mm > 2) {
            jm = mm + 1;
        } else {
            --jy;
            jm = mm + 13;
        }
        jul = (int) (floor(365.25 * jy) + floor(30.6001 * jm) + id + 1720995);
        if (id + 31 * (mm + 12 * iyyy) >= IGREG) { // Test whether to change to
                                                   // Gregorian Cal
            ja = (int) (0.01 * jy); // endar.
            jul += 2 - ja + (int) (0.25 * ja);
        }
        return jul;
    }

    public static void caldat(final int julian, final int[] mm_ref, final int[] id_ref, final int[] iyyy_ref) {
        // Inverse of the function julday given above. Here julian is input
        // as a Julian Day Number, and the routine outputs mm,id, and iyyy as
        // the month, day, and year on which the specified Julian Day
        // started at noon.
        final int IGREG = 2299161;
        int ja, jalpha, jb, jc, jd, je;
        if (julian >= IGREG) { // Cross-over to Gregorian Calendar produces this
                               // correction.
            jalpha = (int) (((double) (julian - 1867216) - 0.25) / 36524.25); 
            ja = julian + 1 + jalpha - (int) (0.25 * jalpha);
        } else if (julian < 0) { // Make day number positive by adding integer
                                 // number of
            // Julian centuries, then subtract them off at the end.
            ja = julian + 36525 * (1 - julian / 36525);
        } else
            ja = julian;
        jb = ja + 1524;
        jc = (int) (6680.0 + ((double) (jb - 2439870) - 122.1) / 365.25);
        jd = (int) (365 * jc + (0.25 * jc));
        je = (int) ((jb - jd) / 30.6001);
        id_ref[0] = jb - jd - (int) (30.6001 * je);
        mm_ref[0] = je - 1;
        if (mm_ref[0] > 12)
            mm_ref[0] -= 12;
        iyyy_ref[0] = jc - 4715;
        if (mm_ref[0] > 2)
            --(iyyy_ref[0]);
        if (iyyy_ref[0] <= 0)
            --(iyyy_ref[0]);
        if (julian < 0)
            iyyy_ref[0] -= 100 * (1 - julian / 36525);
    }

}
