
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Complex.*;

import com.snuggy.nr.chapter17.*;
import com.snuggy.nr.util.*;

public class Hypderiv implements Dtype {

    // Functor to compute derivatives for the hypergeometric equation; see text
    // equation (5.14.4).
    private Complex a, b, c, z0, dz;

    Hypderiv(final Complex aa, final Complex bb, final Complex cc, final Complex z00, final Complex dzz) {
        a = complex(aa);
        b = complex(bb);
        c = complex(cc);
        z0 = complex(z00);
        dz = complex(dzz);
    }

    public void eval(final double s, final double[] yy, final double[] dyyds) {
        Complex z = complex();
        Complex[] y = new Complex[2], dyds = new Complex[2];
        y[0] = complex(yy[0], yy[1]);
        y[1] = complex(yy[2], yy[3]);
        // z=z0+s*dz;
        z = plus(z0, times(s, dz));
        // dyds[0]=y[1]*dz;
        dyds[0] = times(y[1], dz);
        // dyds[1] = (a * b * y[0] - (c - (a + b + 1.0) * z) * y[1]) * dz / (z *
        // (1.0 - z));
        dyds[1] = times(minus(times(a, times(b, y[0])), times(minus(c, times(plus(a, plus(b, 1.0)), z)), y[1])),
                divide(dz, times(z, minus(1.0, z))));
        dyyds[0] = real(dyds[0]);
        dyyds[1] = imag(dyds[0]);
        dyyds[2] = real(dyds[1]);
        dyyds[3] = imag(dyds[1]);
    }

    @Override
    public void jacobian(double x, final double[] y, final double[] dfdx, final double[][] dfdy) throws NRException {
        throw new NRException();
    }

}
