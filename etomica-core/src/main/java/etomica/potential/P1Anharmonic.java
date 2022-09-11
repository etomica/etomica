/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;


public class P1Anharmonic implements IPotential1 {

    private final Space space;
    private double k2 = 1.0, k4 = 1.0;
    private int nBeads;
    private final Vector x0;
    private Tensor eye;

    public P1Anharmonic(Space space, double k2, double k4, int nBeads) {
        this.space = space;
        this.k2 = k2;
        this.k4 = k4;
        this.nBeads = nBeads;
        x0 = space.makeVector();
        eye = space.makeTensor();
        eye.setComponent(0, 0, 1.0);
        eye.setComponent(1, 1, 1.0);
        eye.setComponent(2, 2, 1.0);

    }
    public void setSpringConstants(double k2, double k4) {
        this.k2 = k2;
        this.k4 = k4;
    }
    
    public double[] getSpringConstants() {
        return new double[] {k2, k4};
    }
    
    public void setX0(Vector x0) {
        this.x0.E(x0);
    }
    
    public Vector getX0() {
        return x0;
    }

    public Dimension getX0Dimension() {
        return Length.DIMENSION;
    }

    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION, Length.DIMENSION}, new double[]{1, -2});
    }

    public double u(IAtom atom) {
        Vector dr = space.makeVector();
        dr.Ev1Mv2(atom.getPosition(), x0);
        return 1.0/2.0*k2/nBeads*dr.squared() + 1.0/24.0*k4/nBeads*dr.squared()*dr.squared() ;
    }

    public double udu(IAtom atom, Vector f) {
        Vector dr = space.makeVector();
        dr.Ev1Mv2(atom.getPosition(), x0);
        double u = 1.0/nBeads*(1.0/2.0*k2*dr.squared() + 1.0/24.0*k4*dr.squared()*dr.squared());
        f.PEa1Tv1(-1.0/nBeads*(k2 + k4/6.0*dr.squared()), dr);
        return u;
    }

    @Override
    public Tensor d2u(IAtom atom) {
        Tensor Hii = space.makeTensor();
        Vector dr = space.makeVector();
//        System.out.println(atom + "  " + atom.getPosition());
        dr.Ev1Mv2(atom.getPosition(), x0);
        double dr2 = dr.squared();
        double dr4 = dr2*dr2;
        double du = 1.0/nBeads*(k2*dr2 + 1.0/6.0*k4*dr4);  //rdu: r (k2 r + 1/6 k4 r^3)
        double d2u = 1.0/nBeads*(k2*dr2 + 1.0/2.0*k4*dr4); //r2d2u: r^2 (k2 + 1/2 k4 r^2)
        Hii.Ev1v2(dr, dr);
        Hii.TE((d2u-du)/dr4);
        Hii.PEa1Tt1(du/dr2, eye);
        return Hii;
    }

}
   
