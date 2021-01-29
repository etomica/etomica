/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Exponential-6 atom-atom repulsion-dispersion potential. Given formula:
 * <p>
 * U(r) = A*exp(-r/B) - C /r^6
 *
 * @author Tai Tan
 */

public class P2Exp6 extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, double AA, double BB, double CC, TruncationFactory tf) {
        return tf.make(new P2Exp6(space, AA, BB, CC));
    }

    public P2Exp6(Space _space) {
        // these defaults probably aren't appropriate -- need to develop A,B,C
        // from default size, well depth, and well extent (which doesn't exist!
        // maybe potl cutoff?)
        this(_space, 1.0, 1.0, 1.0);

    }

    public P2Exp6(Space _space, double AA, double BB, double CC) {
        super(_space);
        dr01 = space.makeVector();
        setA(AA);
        setB(BB);
        setC(CC);
    }

    /**
     * The energy u
     */
    public double u(double r2) {
        double r = Math.sqrt(r2);

        if (r < 3 * BB) {
            return Double.POSITIVE_INFINITY;
        }
        return AA * Math.exp(-r / BB) - CC / (r2 * r2 * r2);
        
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        return -(AA / BB) * r * Math.exp(-r / BB) + 6 * CC / (r2 * r2 * r2);

    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation: r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        return (AA / (BB * BB)) * r2 * Math.exp(-r / BB) - 42 * CC
                / (r2 * r2 * r2);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) { // need long range correction!!!!
        throw new MethodNotImplementedException("Integral for long-range correction for Exp-6 not yet implemented");
    }

    public double getA() {
        return AA;
    }

    public final void setA(double a) {
        AA = a;
    }

    public double getB() {
        return BB;
    }

    public final void setB(double b) {
        BB = b;
    }

    public double getC() {
        return CC;
    }

    public final void setC(double c) {
        CC = c;
    }

    public Dimension getADimension() {
        return Energy.DIMENSION;
    }

    public Dimension getBDimension() {
        return Length.DIMENSION;
    }

    public Dimension getCDimension() {
        return new CompoundDimension(new Dimension[] {Energy.DIMENSION, Length.DIMENSION}, new double[] {1.0, 6.0});
    }

    private double AA, BB, CC;
    protected final Vector dr01;

    private static final long serialVersionUID = 1L;
}
