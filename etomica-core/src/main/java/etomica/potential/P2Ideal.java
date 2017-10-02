/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Ideal-gas two-body potential, which defines no interactions and zero energy
 * for all pairs given to it.
 * <p> 
 * Useful as a placeholder where a potential is expected but it is desired to 
 * not have the atoms interact.
 */
public class P2Ideal extends Potential2 implements Potential2Soft,
        Potential2Spherical, PotentialHard {

    public P2Ideal(Space space) {
        super(space);
        zeroVector = new Vector[1];
        zeroVector[0] = space.makeVector();
        zeroTensor = space.makeTensor();
    }
    /**
     * Does nothing.
     */
    public void setBox(Box box) {
    }

    public void setRange(double range){
        this.range = range;
    }

    /**
     * Returns zero.
     */
    public double getRange() {
        return range;
    }

    /**
     * Returns zero.
     */
    public double energy(IAtomList atoms) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double hyperVirial(IAtomList pair) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double virial(IAtomList pair) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double integral(double rC) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double u(double r2) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double du(double r2) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double lastCollisionVirial() {
        return 0;
    }

    /**
     * Returns a tensor of zeros.
     */
    public Tensor lastCollisionVirialTensor() {
        zeroTensor.E(0.0);
        return zeroTensor;
    }

    /**
     * Does nothing.
     */
    public void bump(IAtomList atom, double falseTime) {
    }

    /**
     * Returns Double.POSITIVE_INFINITY.
     */
    public double collisionTime(IAtomList atom, double falseTime) {
        return Double.POSITIVE_INFINITY;
    }

    /**
     * Returns zero.
     */
    public double energyChange() {
        return 0;
    }

    /**
     * Returns a zero vector.
     */
    public Vector[] gradient(IAtomList atoms) {
        zeroVector[0].E(0.0);
        return zeroVector;
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
        

    private static final long serialVersionUID = 1L;
    private final Vector[] zeroVector;
    private final Tensor zeroTensor;
    protected double range;
}
