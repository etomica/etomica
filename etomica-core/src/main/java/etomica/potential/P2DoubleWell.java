/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.Debug;

/**
 * Double-well potential.
 * Energy is -epsilonCore if distance is less than sigma,
 * is -epsilonWell if less than lambda*sigma and not overlapping,
 * and is zero otherwise.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 */
public class P2DoubleWell extends Potential2HardSpherical {

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilonCore, epsilonWell;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected Vector dv;

    public P2DoubleWell(Space space) {
        this(space, 1.0, 2.0, Double.POSITIVE_INFINITY, 1.0);
    }

    public P2DoubleWell(Space space, double coreDiameter, double lambda, double epsilonCore, double epsilonWell) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilonCore(epsilonCore);
        setEpsilonWell(epsilonWell);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
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
        double ke = bij*bij*reduced_m/(2*r2);
        if(2*r2 < (wellDiameterSquared+coreDiameterSquared)) {   // Hard-core collision
            if (Debug.ON && Math.abs(r2 - coreDiameterSquared)/coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+coreDiameterSquared);
            }
            if (bij > 0) {
                if (ke < (epsilonCore-epsilonWell)) { // Not enough kinetic energy to escape inner well
                    lastCollisionVirial = 2*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {  // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*(epsilonCore-epsilonWell)/reduced_m));
                    nudge = eps;
                    lastEnergyChange = (epsilonCore-epsilonWell);
                }
            }
            else if (ke > -(epsilonCore-epsilonWell)) { // approach / capture
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*(epsilonCore-epsilonWell)/reduced_m));
                nudge = -eps;
                lastEnergyChange = -epsilonCore;
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
                if(ke < epsilonWell) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilonWell/reduced_m));
                    nudge = eps;
                    lastEnergyChange = epsilonWell;
                }
            }
            else if(ke > -epsilonWell) {   // Approach/capture
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilonWell/reduced_m));
                nudge = -eps;
                lastEnergyChange = -epsilonWell;
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
        if (rm0 > 0) {
            atom0.getPosition().PEa1Tv1(-nudge,dr);
        }
        if (rm1 > 0) {
            atom1.getPosition().PEa1Tv1(nudge,dr);
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
        if (coord0.getParentGroup().getType().getIndex()+coord0.getParentGroup().getType().getIndex()>0) {
            System.out.println("here we go");
        }
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double time = Double.POSITIVE_INFINITY;

        if(r2 < wellDiameterSquared) {  // Already inside wells
            if (r2 < coreDiameterSquared) { // inside core
                double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                time = (-bij + Math.sqrt(discr))/v2;
            }
            else {
                if(bij < 0.0) {    // Check for hard-core collision
    
                    double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                    if(discr > 0) {  // Hard cores collide next
                        time = (-bij - Math.sqrt(discr))/v2;
                        if (time<0) {
                            throw new RuntimeException("oops");
                        }
                    }
                    else {           // Moving toward each other, but wells collide next
                        discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                        time = (-bij + Math.sqrt(discr))/v2;
                        if (time<0) {
                            throw new RuntimeException("oops");
                        }
                    }
                }
                else {           // Moving away from each other, wells collide next
                    double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                    time = (-bij + Math.sqrt(discr))/v2;
                }
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
        if (time < 0) {
            throw new RuntimeException("oops");
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
        if (r2 > coreDiameterSquared) return -epsilonWell;
        return -epsilonCore;
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
    public double getEpsilonCore() {return epsilonCore;}
    public double getEpsilonWell() {return epsilonWell;}
    /**
     * Accessor method for depth of well
     */
    public void setEpsilonCore(double epsCore) {
        epsilonCore = epsCore;
    }
    public void setEpsilonWell(double epsWell) {
        epsilonWell = epsWell;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

}
  
