package etomica.potential;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Ideal-gas two-body potential, which defines no interactions and zero energy
 * for all pairs given to it.
 * <p> 
 * Useful as a placeholder where a potential is expected but it is desired to 
 * not have the atoms interact.
 */

/*
 * History
 * Created on Jul 21, 2005 by kofke
 */
public class P2Ideal extends Potential2 implements Potential2Soft,
        Potential2Spherical, PotentialHard, PotentialSoft {

    public P2Ideal(Space space) {
        super(space);
        zeroVector = new IVector[1];
        zeroVector[0] = space.makeVector();
        zeroTensor = space.makeTensor();
    }
    /**
     * Does nothing.
     */
    public void setPhase(Phase phase) {
    }

    /**
     * Returns zero.
     */
    public double getRange() {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double energy(AtomSet atoms) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double hyperVirial(AtomSet pair) {
        return 0;
    }

    /**
     * Returns zero.
     */
    public double virial(AtomSet pair) {
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
    public void bump(AtomSet atom, double falseTime) {
    }

    /**
     * Returns Double.POSITIVE_INFINITY.
     */
    public double collisionTime(AtomSet atom, double falseTime) {
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
    public IVector[] gradient(AtomSet atoms) {
        zeroVector[0].E(0.0);
        return zeroVector;
    }

    private final IVector[] zeroVector;
    private final Tensor zeroTensor;
}
