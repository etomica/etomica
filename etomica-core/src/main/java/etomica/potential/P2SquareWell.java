/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.Debug;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 * 
 * This implementation determines which atoms are interacting (in the well) by
 * checking their distance.  This is prone to numerical roundoff issues, and so
 * atoms are bumped off the well after capture and escape.  This can cause
 * issues related to 3-particle collisions.  P2SquareWellRobust can be used
 * instead to avoid such issues.
 */
public class P2SquareWell extends Potential2HardSpherical {

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected Vector dv;
    protected final boolean ignoreOverlap;

    public P2SquareWell(Space space) {
        this(space, 1.0, 2.0, 1.0, false);
    }

    public P2SquareWell(Space space, double coreDiameter, double lambda, double epsilon, boolean ignoreOverlap) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
        this.ignoreOverlap = ignoreOverlap;
    }

    public double getRange() {
        return wellDiameter;
    }
    
    /**
     * Implements collision dynamics between two square-well atoms.
     * Includes all possibilities involving collision of hard cores, and collision of wells
     * both approaching and diverging
     */
    public void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double eps = 1.0e-10;
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reduced_m = 1.0/(rm0+rm1);
        double nudge = 0;
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared) || lambda == 1) {   // Hard-core collision
            if (Debug.ON && !ignoreOverlap && Math.abs(r2 - coreDiameterSquared)/coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+coreDiameterSquared);
            }
            if (bij > 0) {
                lastCollisionVirial = 0;
                nudge = eps;
            }
            else {
                lastCollisionVirial = 2.0*reduced_m*bij;
            }
            lastEnergyChange = 0.0;
        }
        else {    // Well collision
            if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+wellDiameterSquared);
            }
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(bij > 0.0) {         // Separating
                if(ke < epsilon) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
                    nudge = eps;
                    lastEnergyChange = epsilon;
                }
            }
            else if(ke > -epsilon) {   // Approach/capture
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
                nudge = -eps;
                lastEnergyChange = -epsilon;
            }
            else {                     // Not enough kinetic energy to overcome square-shoulder
                lastCollisionVirial = 2.0*reduced_m*bij;
                nudge = eps;
                lastEnergyChange = 0.0;
            }
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1( rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);
        if(nudge != 0) {
            if (rm0 > 0) {
                atom0.getPosition().PEa1Tv1(-nudge,dr);
            }
            if (rm1 > 0) {
                atom1.getPosition().PEa1Tv1(nudge,dr);
            }
        }
    }//end of bump method

    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic coord0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic coord1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double time = Double.POSITIVE_INFINITY;

        if(r2 < wellDiameterSquared) {  // Already inside wells

            if(bij < 0.0) {    // Check for hard-core collision
                if(ignoreOverlap && r2 < coreDiameterSquared) {   // Inside core; collision now
                    return falseTime+0.001*Math.sqrt(dr.squared())/Math.sqrt(v2);
                }

                double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                if(discr > 0) {  // Hard cores collide next
                    time = (-bij - Math.sqrt(discr))/v2;
                }
                else {           // Moving toward each other, but wells collide next
                    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                    time = (-bij + Math.sqrt(discr))/v2;
                }
            }
            else {           // Moving away from each other, wells collide next
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                time = (-bij + Math.sqrt(discr))/v2;
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                if(discr > 0) {
                    time = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.allAtoms(pair)) || time < 0)) {
            System.out.println(pair+" r2 "+r2+" bij "+bij+" time "+(time+falseTime));
        }
        return time + falseTime;
    }

  /**
   * Returns infinity if overlapping, -epsilon if otherwise less than well diameter, or zero if neither.
   */
    public double u(double r2) {
        if (r2 > wellDiameterSquared) return 0.0;
        if (r2 > coreDiameterSquared) return -epsilon;
        return Double.POSITIVE_INFINITY;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return coreDiameter;}
    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        if (coreDiameter==0) wellDiameter = lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getCoreDiameterDimension() {return Length.DIMENSION;}

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {return lambda;}
    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when 
     * this is changed
     */
    public void setLambda(double lam) {
        if (lam < 1.0) throw new IllegalArgumentException("Square-well lambda must be greater than 1.0");
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        if (coreDiameter==0) wellDiameter = lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getLambdaDimension() {return Null.DIMENSION;}

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for depth of well
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

}
  
