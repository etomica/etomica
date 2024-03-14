/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;


public class P1AnharmonicTIA implements IPotential1 {

    private final Space space;
    private double k2 = 1.0, k4 = 1.0;
    private final Vector x0;
    private int nBeads;
    private double facTIA;
    private double mOmegaN2;


    public P1AnharmonicTIA(Space space, double k2, double k4, int nBeads, double mOmegaN2, double facTIA) {
        this.space = space;
        this.k2 = k2;
        this.k4 = k4;
        this.nBeads = nBeads;
        this.facTIA = facTIA;
        this.mOmegaN2 = mOmegaN2;
        x0 = space.makeVector();
    }
    public void setSpringConstants(double k2, double k4) {
        this.k2 = k2;
        this.k4 = k2;
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
        double dr2 = dr.squared();
        double U =  1.0/2.0*k2/nBeads*dr2 + 1.0/24.0*k4/nBeads*dr2*dr2
                + facTIA*1.0/24.0/(mOmegaN2*nBeads)/nBeads*(k2*k2*dr2 + 1.0/3.0*k2*k4*dr2*dr2 + 1.0/36.0*k4*k4*dr2*dr2*dr2);
        return U;
    }

    public double udu(IAtom atom, Vector f) {
        Vector dr = space.makeVector();
        dr.Ev1Mv2(atom.getPosition(), x0);
        double dr2 = dr.squared();
        double U =  1.0/2.0*k2/nBeads*dr2 + 1.0/24.0*k4/nBeads*dr2*dr2
                + facTIA*1.0/24.0/(mOmegaN2*nBeads)/nBeads*(k2*k2*dr2 + 1.0/3.0*k2*k4*dr2*dr2 + 1.0/36.0*k4*k4*dr2*dr2*dr2);

        f.PEa1Tv1(-k2/nBeads - k4/nBeads/6.0*dr2, dr);
        f.PEa1Tv1(-facTIA/24.0/(mOmegaN2*nBeads)/nBeads*(2.0*k2*k2 + 4.0/3.0*k2*k4*dr2 + 1.0/6.0*k4*k4*dr2*dr2), dr);
        return U;
    }

}
   
