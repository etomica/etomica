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
import etomica.util.Debug;
/**
 * Basic penetrable-sphere potential.
 * Energy is epsilon if spheres overlap, and is zero otherwise.  
 * Core diameter describes height of penetrable core.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce a square-well potential with no hard core. 
 */
public class P2PenetrableSphere extends Potential2HardSpherical {

    private static final long serialVersionUID = 1L;
    protected double coreDiameter, coreDiameterSquared;
    protected double epsilon;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected Vector dv;

    /**
     * Constructor with default values of unity for core diameter and core energy.
     */
    public P2PenetrableSphere(Space space) {
        this(space, 1.0, 1.0);
    }

    public P2PenetrableSphere(Space space, double coreDiameter, double epsilon) {
        super(space);
        setCoreDiameter(coreDiameter);
        setEpsilon(epsilon);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
    }

    public double getRange() {
        return coreDiameter;
    }
    
    /**
     * Implements collision dynamics between two penetrable-sphere atoms.
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
        if(bij > 0.0) {         // Separating
            if(ke < -epsilon) {     // Core is a well -- Not enough kinetic energy to escape//
                lastCollisionVirial = 2.0*reduced_m*bij;
                nudge = -eps;
                lastEnergyChange = 0.0;
            }
            else {                 // Core is a shoulder, or is a well and KE is enough to escape -- separate
                lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij + 2.0*r2*epsilon/reduced_m));//
                nudge = eps;
                lastEnergyChange = -epsilon;
            }
        }
        else if(ke > epsilon) {   // Approach and KE is enough to enter core
            lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij-2.0*r2*epsilon/reduced_m));
            nudge = -eps;
            lastEnergyChange = epsilon;
        }
        else {                     // Not enough kinetic energy to overcome square-shoulder
            lastCollisionVirial = 2.0*reduced_m*bij;
            nudge = eps;
            lastEnergyChange = 0.0;
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
        if(r2 < coreDiameterSquared) {  // Already inside PS; MUST leave PS with all bij and discr!
            double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
            time = (-bij + Math.sqrt(discr))/v2;
        } else if (bij<0.0) {// Outside PS ; look for collision at PS
            double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
            if (discr > 0.0 ){          
                    time = (-bij - Math.sqrt(discr))/v2;
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.allAtoms(pair)) || time < 0)) {
            System.out.println(pair+" r2 "+r2+" bij "+bij+" time "+(time+falseTime));
        }
        return time + falseTime;
    }

  /**
   * Returns epsilon if overlapping, or zero if not.
   */
    public double u(double r2) {
        if (r2 > coreDiameterSquared) return 0;
        return epsilon;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return coreDiameter;}
    
    /**
     * Mutator method for core diameter.
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        coreDiameter = c;
        coreDiameterSquared = c*c;

    }
    public Dimension getCoreDiameterDimension() {return Length.DIMENSION;}

    /**
     * Mutator method for height of core.
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for height of core. Positive value corresponds to a positive core height.
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

}
  
