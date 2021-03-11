/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * pseudo-potential for a "collision" time to update colliders for periodic boundaries
 */
 
public class P1HardPeriodic extends Potential1 implements PotentialHard {

    /**
     * Returns an instance of P1HardPeriodic with sigma = NaN.  call setSigma
     * to set the value you want.
     */
    public P1HardPeriodic(Space space) {
        this(space, Double.NaN);
        // use NaN so they'll have to call setSigma later
    }

    /**
     * Returns an instance of P1HardPeriodic with the given value of sigma (the
     * maximum distance between two atoms where they interact)
     */
    public P1HardPeriodic(Space space, double sigma) {
        super(space);
        this.sigma = sigma;
    }
    
    /**
     * Returns zero.
     */
    public double energy(IAtomList a) {
        return 0.0;
    }
     
    /**
     * Returns zero.
     */
    public double energyChange() {
        return 0.0;
    }
    
    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.get(0);
        Vector v = atom.getVelocity();
        Vector dim = boundary.getBoxSize();
        double tmin = Double.POSITIVE_INFINITY;
        double d2 = 2.0*sigma;
        int D = dim.getD();
        // 4*(L/4 - sigma/2) = L - 2sigma
        for(int i=0; i<D; i++) {
            double t = (dim.getX(i)-d2)/v.getX(i);
            t = (t < 0) ? -t : t;//abs
            tmin = (t < tmin) ? t : tmin;
        }
        return 0.25*tmin + falseTime;
    }
    
    public void setSigma(double newSigma) {
        sigma = newSigma;
    }
    
    public double getSgima() {
        return sigma;
    }
    
    /**
     * Performs no action.
     */
    public void bump(IAtomList a, double falseTime) { }
    
    /**
     * Returns zero.
     */
    public double lastCollisionVirial() {return 0;}
    
    /**
     * Returns null.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
    private static final long serialVersionUID = 1L;
    protected double sigma;
}
   
