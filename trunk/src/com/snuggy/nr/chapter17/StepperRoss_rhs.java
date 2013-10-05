
package com.snuggy.nr.chapter17;

public class StepperRoss_rhs implements Dtype {
    @Override
    public void eval(final double x,final double[] y,final double[] dydx) {
      dydx[0]= -0.013*y[0]-1000.0*y[0]*y[2];
      dydx[1]= -2500.0*y[1]*y[2];
      dydx[2]= -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
    }
    @Override
    public void jacobian(final double x,final double[] y,final double[] dfdx,final double[][] dfdy) {
      int n=y.length;
      for (int i=0;i<n;i++) dfdx[i]=0.0;
      dfdy[0][0]= -0.013-1000.0*y[2];
      dfdy[0][1]= 0.0;
      dfdy[0][2]= -1000.0*y[0];
      dfdy[1][0]= 0.0;
      dfdy[1][1]= -2500.0*y[2];
      dfdy[1][2]= -2500.0*y[1];
      dfdy[2][0]= -0.013-1000.0*y[2];
      dfdy[2][1]= -2500.0*y[2];
      dfdy[2][2]= -1000.0*y[0]-2500.0*y[1];
    }
  }