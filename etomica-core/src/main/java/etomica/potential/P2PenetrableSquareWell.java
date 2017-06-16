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
 * Square-well potential with a penetrable core.
 * Energy is epsilonCore if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of penetrable core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 */
public class P2PenetrableSquareWell extends Potential2HardSpherical {

    private static final long serialVersionUID = 1L;
    protected double coreDiameterSquared;
    protected double wellDiameterSquared;
    protected double lambdaSquared; //wellDiameter = coreDiameter * lambda
    protected double epsilon, epsilonCore, deltaEps;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected Vector dv;
    protected final boolean ignoreOverlap;

    public P2PenetrableSquareWell(Space space) {
        this(space, 1.0, 2.0, 1.0, false);
    }

    public P2PenetrableSquareWell(Space space, double coreDiameter, double lambda, double epsilon, boolean ignoreOverlap) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
        this.ignoreOverlap = ignoreOverlap;
    }

    public double getRange() {
        return Math.sqrt(wellDiameterSquared);
    }
    
    /**
     * Implements collision dynamics between two square-well atoms.
     * Includes all possibilities involving collision of cores, and collision of wells
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
        double ke = bij*bij*reduced_m/(2.0*r2);
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared) || lambdaSquared == 1) {   // Hard-core collision
            if (Debug.ON && !ignoreOverlap && Math.abs(r2 - coreDiameterSquared)/coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+coreDiameterSquared);
            }
            if(bij > 0.0) {         // Separating
                if(ke < -deltaEps) {     // Core is a well -- Not enough kinetic energy to escape//
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Core is a shoulder, or is a well and KE is enough to escape -- separate
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij + 2.0*r2*deltaEps/reduced_m));//
                    nudge = eps;
                    lastEnergyChange = -deltaEps;
                }
            }
            else if(ke > deltaEps) {   // Approach and KE is enough to enter core
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij-2.0*r2*deltaEps/reduced_m));
                nudge = -eps;
                lastEnergyChange = deltaEps;
            }
            else {                     // Not enough kinetic energy to overcome square-shoulder
                lastCollisionVirial = 2.0*reduced_m*bij;
                nudge = eps;
                lastEnergyChange = 0.0;
            }
        }
        else {    // Well collision
            if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+wellDiameterSquared);
            }
            // ke is kinetic energy due to components of velocity
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
        
        if(r2 < coreDiameterSquared) {  // Already inside penetrable core
            double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
            time = (-bij + Math.sqrt(discr))/v2;
            
        } else if(r2 < wellDiameterSquared) {  // Already inside wells           
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
        } else {              // Outside wells; look for collision at well
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
        return epsilonCore;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return Math.sqrt(coreDiameterSquared);}
    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        setCoreDiameterSquared(c*c);
    }
    public Dimension getCoreDiameterDimension() {return Length.DIMENSION;}
    
    public void setCoreDiameterSquared(double c2) {
        coreDiameterSquared = c2;
        wellDiameterSquared = c2*lambdaSquared;
    }

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {return Math.sqrt(lambdaSquared);}
    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when 
     * this is changed
     */
    public void setLambda(double lam) {
        if (lam < 1.0) throw new IllegalArgumentException("Square-well lambda must be greater than 1.0");
        lambdaSquared = lam*lam;
        wellDiameterSquared = coreDiameterSquared*lambdaSquared;
    }
    public Dimension getLambdaDimension() {return Null.DIMENSION;}

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {return epsilon;}
    
    /**
     * Mutator method for depth of well. Positive value corresponds to depth of the well, and a negative well energy.
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
        deltaEps = epsilonCore + epsilon;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

    public double getEpsilonCore() {
        return epsilonCore;
    }

    /**
     * Sets height of penetrable core.  Positive value corresponds to a positive height of the core energy.
     */
    public void setEpsilonCore(double epsilsonCore) {
        this.epsilonCore = epsilsonCore;
        deltaEps = epsilonCore + epsilon;
    }
    public Dimension getEpsilonCoreDimension() {return Energy.DIMENSION;}

}
  
