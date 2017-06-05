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
 * Potential with a well having a finite potential on either side. Similar to P2HardBond,
 * but allows for separations greater or less than the well separation. 
 */
public class P2PenetrableBond extends Potential2HardSpherical {

    public P2PenetrableBond(Space space) {
        this(space, 1.0, 0.15, false);
    }

    public P2PenetrableBond(Space space, double bondLength, double bondDelta, boolean ignoreOverlap) {
        super(space);
        setBondLength(bondLength);
        setBondDelta(bondDelta);
        lastCollisionVirialTensor = space.makeTensor();
        dv = space.makeVector();
        this.ignoreOverlap = ignoreOverlap;
    }

    /**
     * Accessor method for the bond length
     */
    public final double getBondLength() {
        return bondLength;
    }

    /**
     * Accessor method for the bond extension factor
     */
    public final double getBondDelta() {
        return bondDelta;
    }

    /**
     * Setter method for the bond length
     */
    public final void setBondLength(double l) {
        bondLength = l;
        double minBondLength = bondLength;
        maxBondLength = bondLength + bondDelta;
        minBondLengthSquared = minBondLength * minBondLength;
        maxBondLengthSquared = maxBondLength * maxBondLength;
    }

    /**
     * Sets the bond extension factor (max = length * (1+delta))
     */
    public final void setBondDelta(double d) {
        if(d == 0) throw new IllegalArgumentException("P2PenetrableBond is not set to handle bondDelta = 0");
        bondDelta = d;
        setBondLength(bondLength);
    }

    public final Dimension getBondLengthDimension() {
        return Length.DIMENSION;
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
        if(2*r2 < (minBondLengthSquared+maxBondLengthSquared)) {   // Hard-core collision
            if (Debug.ON && !ignoreOverlap && Math.abs(r2 - minBondLengthSquared)/minBondLengthSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+minBondLengthSquared);
            }
            
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(bij < 0.0) {         // Approaching
                if(ke < epsInner) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij + Math.sqrt(bij*bij - 2.0*r2*epsInner/reduced_m));
                    nudge = -eps;
                    lastEnergyChange = epsInner;
                }
            }
                        
            else if(ke > -epsInner) {   // Separation
                double epsEffective = epsInner;
                if(Double.isInfinite(epsInner)) epsEffective = 0.0;//ignore effect of capture if from infinite potential
                lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij+2.0*r2*epsEffective/reduced_m));
                nudge = eps;
                lastEnergyChange = -epsEffective;
            }
            else {                     // Not enough kinetic energy to overcome square-shoulder
                lastCollisionVirial = 2.0*reduced_m*bij;
                nudge = -eps;
                lastEnergyChange = 0.0;
            }

        }
        else {    // Well collision
            if (Debug.ON && Math.abs(r2 - maxBondLengthSquared)/maxBondLengthSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+maxBondLengthSquared);
            }
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(bij > 0.0) {         // Separating
                if(ke < epsOuter) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsOuter/reduced_m));
                    nudge = eps;
                    lastEnergyChange = epsOuter;
                }
            }
            else if(ke > -epsOuter) {   // Approach/capture
                double epsEffective = epsOuter;
                if(Double.isInfinite(epsOuter)) epsEffective = 0.0;//ignore effect of capture if from infinite potential
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsEffective/reduced_m));
                nudge = -eps;
                lastEnergyChange = -epsEffective;
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


    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

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
        
        if(r2 < maxBondLengthSquared) {  // Already inside wells

            if(bij < 0.0) {    // inside wells, moving toward each other
                if(r2 < minBondLengthSquared) {   //inside core, moving toward each other
                    double discr = bij*bij - v2 * ( r2 - minBondLengthSquared );
                    time = (-bij + Math.sqrt(discr))/v2;
                    
//                    return falseTime+0.001*Math.sqrt(dr.squared())/Math.sqrt(v2);
                } else { //inside well but not inside core, moving toward each other
                    double discr = bij*bij - v2 * ( r2 - minBondLengthSquared );
                    if(discr > 0) {  // Hard cores collide next
                        time = (-bij - Math.sqrt(discr))/v2;
                    }
                    else {           // Moving toward each other, but wells collide next
                        discr = bij*bij - v2 * ( r2 - maxBondLengthSquared );
                        time = (-bij + Math.sqrt(discr))/v2;
                    }
                }
            }  
            //moving apart
            else {           
                double discr = 0.0;
                if(r2 < minBondLengthSquared) {//inside core, moving away from each other
                    discr = bij*bij - v2 * ( r2 - minBondLengthSquared );  // This is always > 0
                } else {    //inside well, moving away
                    discr = bij*bij - v2 * ( r2 - maxBondLengthSquared );  // This is always > 0
                }
                time = (-bij + Math.sqrt(discr))/v2;
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                double discr = bij*bij - v2 * ( r2 - maxBondLengthSquared );
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
     * Returns 0 if the bond is within the required distance, infinity if not.
     */
    public double u(double r2) {
        if (r2 < minBondLengthSquared) return epsInner;
        if (r2 < maxBondLengthSquared) return 0.0;
        return epsOuter;
    }
    
    /**
     * Returns infinity.  
     */
    public double getRange() {
        return bondLength + bondDelta;
    }
    
    
    
    public double getEpsInner() {
        return epsInner;
    }

    public void setEpsInner(double epsInner) {
        this.epsInner = epsInner;
    }

    public double getEpsOuter() {
        return epsOuter;
    }

    public void setEpsOuter(double epsOuter) {
        this.epsOuter = epsOuter;
    }
    
    public Dimension getEpsInnerDimension() {
        return Energy.DIMENSION;
    }

    public Dimension getEpsOuterDimension() {
        return Energy.DIMENSION;
    }

    public double energyChange() {return lastEnergyChange;}

    private static final long serialVersionUID = 1L;
    private double minBondLengthSquared;
    private double lastEnergyChange = 0.0;
    private double maxBondLength, maxBondLengthSquared;
    private double bondLength;
    private double bondDelta;
    private double epsInner, epsOuter;
    private double lastCollisionVirial = 0.0;
    private double lastCollisionVirialr2 = 0.0;
    private boolean ignoreOverlap;
    private final Vector dv;
    private final Tensor lastCollisionVirialTensor;
}
