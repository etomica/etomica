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


public class P1Anharmonic234 implements IPotential1 {

    private final Space space;
    private double k2 = 1.0, k4 = 1.0, k3 = 1.0;
    private final Vector x0;
    private Tensor eye;

    public P1Anharmonic234(Space space, double k2, double k3, double k4) {
        this.space = space;
        this.k2 = k2;
        this.k3 = k3;
        this.k4 = k4;
        x0 = space.makeVector();
        eye = space.makeTensor();
        eye.setComponent(0, 0, 1.0);
        eye.setComponent(1, 1, 1.0);
        eye.setComponent(2, 2, 1.0);

    }
    public void setSpringConstants(double k2, double k3, double k4) {
        this.k2 = k2;
        this.k3 = k3;
        this.k4 = k4;
    }
    
    public double[] getSpringConstants() {
        return new double[] {k2, k3, k4};
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
        Vector dx = space.makeVector();
        dx.Ev1Mv2(atom.getPosition(), x0);
        double dx2 = dx.getX(0)*dx.getX(0);
        double dx3 = dx2*dx.getX(0);
        double dx4 = dx2*dx2;
        return 1.0/2.0*k2*dx2 + k3*dx3 + k4*dx4;
    }

    public double udu(IAtom atom, Vector f) {
        Vector dx = space.makeVector();
        dx.Ev1Mv2(atom.getPosition(), x0);
        double dx2 = dx.getX(0)*dx.getX(0);
        double dx3 = dx2*dx.getX(0);
        double dx4 = dx2*dx2;
        double u = 1.0/2.0*k2*dx2 + k3*dx3 + k4*dx4;
        f.PEa1Tv1(-k2 - 3*k3*dx.getX(0) - 4*k4*dx2, dx);
        return u;
    }

    @Override
    public Tensor d2u(IAtom atom) {
        Tensor Hii = space.makeTensor();
        Vector dx = space.makeVector();
        dx.Ev1Mv2(atom.getPosition(), x0);
        double dx2 = dx.getX(0)*dx.getX(0);
        double dx3 = dx2*dx.getX(0);
        double dx4 = dx2*dx2;
        double du = k2*dx2 + 3*k3*dx3 + 4*k4*dx4;//xdu
        double d2u = k2*dx2 + 6*k3*dx3 + 12*k4*dx4;
        Hii.Ev1v2(dx, dx);
        Hii.TE((d2u-du)/dx4);
        Hii.PEa1Tt1(du/dx2, eye);
        return Hii;
    }

}
   
