package etomica.potential;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public abstract class Potential2SoftSpherical extends Potential2 implements Potential2Soft, Potential2Spherical {
   
    public Potential2SoftSpherical(Space space) {
        super(space);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
    }
        
    /**
     * The pair energy u(r^2).
     * @param the square of the distance between the particles.
     */
    public abstract double u(double r2);
        
    /**
     * The derivative of the pair energy, times the separation r: r du/dr.
     */
    public abstract double du(double r2);
        
    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public abstract double d2u(double r2);
        
    /**
     * Integral of the potential, used to evaluate corrections for potential truncation.
     * Specifically, this is the integral from rC (the argument) to infinity of
     * u(r) A r^(D-1), where D is the spatial dimension, and A is the area of a unit
     * sphere in D dimensions.  Normally, the long-range potential correction would be obtained
     * by multiplying this quantity by the pair density nPairs/V, where nPairs is the number of pairs of atoms
     * affected by this potential, and V is the volume they occupy.
     */
    public abstract double uInt(double rC);
    
    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
        nearestImageTransformer.nearestImage(dr);
        return u(dr.squared());
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
        nearestImageTransformer.nearestImage(dr);
        return du(dr.squared());
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        return d2u(r2) + du(r2);
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector[] gradient(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        gradient[1].Ea1Tv1(du(r2)/r2,dr);
        gradient[0].Ea1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    /**
     * Same as uInt.
     */
    public double integral(double rC) {
        return uInt(rC);
    }
    
    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setPhase(Phase phase) {
        nearestImageTransformer = phase.getBoundary();
    }

    private final Vector[] gradient;
    protected NearestImageTransformer nearestImageTransformer;
    protected final Vector dr;
    
}//end of Potential2SoftSpherical
