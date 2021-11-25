/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Simple electrostatic potential class.
 *
 * @author Andrew Schultz
 */
public class P2Electrostatic implements Potential2Soft {

    public static Potential2Soft makeTruncated(double q1, double q2, TruncationFactory tf) {
        return tf.make(new P2Electrostatic(q1, q2));
    }

    public P2Electrostatic(double q1, double q2) {
        setCharge1(q1);
        setCharge2(q2);
    }

    public double d2u(double r2) {
        return +2 * u(r2);
    }

    public double du(double r2) {
        return -u(r2);
    }

    public double u(double r2) {
        return qq/Math.sqrt(r2);
    }

    /**
     * Returns the charge on the first type of atom.
     */
    public double getCharge1() {
        return q1;
    }
    
    /**
     * Sets the charge on the first type of atom.
     */
    public void setCharge1(double newCharge1) {
        q1 = newCharge1;
        qq = q1*q2;
    }

    /**
     * Returns the charge on the second type of atom.
     */
    public double getCharge2() {
        return q2;
    }
    
    /**
     * Sets the charge on the second type of atom.
     */
    public void setCharge2(double newCharge2) {
        q2 = newCharge2;
        qq = q1*q2;
    }
    
    private static final long serialVersionUID = 1L;
    protected double q1, q2, qq;
}
