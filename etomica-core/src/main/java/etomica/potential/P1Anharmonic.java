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

/**
 * Potential in which attaches a harmonic spring between each affected atom and
 * the nearest boundary in each direction.
 * <p>
 * This class has not been used or checked for correctness.
 *
 * @author David Kofke
 */

public class P1Anharmonic implements IPotential1 {

    private final Space space;
    private double k2 = 1.0, k4 = 1.0;
    private final Vector x0;

    public P1Anharmonic(Space space, double springConstant1, double springConstant2) {
        this.space = space;
        this.k2 = springConstant1;
        this.k4 = springConstant2;
        x0 = space.makeVector();
    }
    public void setSpringConstants(double springConstant1, double springConstant2) {
        k2 = springConstant1;
        k4 = springConstant2;
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
        return 1.0/2.0*k2*atom.getPosition().Mv1Squared(x0) +
               1.0/24.0*k4*atom.getPosition().Mv1Squared(x0)*atom.getPosition().Mv1Squared(x0) ;
    }

    public double udu(IAtom atom, Vector f) {
        Vector dr = space.makeVector();
        dr.Ev1Mv2(atom.getPosition(), x0);
        double u = 1.0/2.0*k2*dr.squared()+1.0/24.0*k4*dr.squared()*dr.squared();
        f.PEa1Tv1(-k2 - k4/6.0*dr.squared(), dr);
        return u;
    }

}
   
