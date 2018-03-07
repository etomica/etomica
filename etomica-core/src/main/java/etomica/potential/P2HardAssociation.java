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

/**
 * Purely attractive square-well potential with no repulsive core.  Similar
 * to P2Tether, but permits atoms to separate.
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
public class P2HardAssociation extends Potential2HardSpherical {

    public P2HardAssociation(Space space) {
        this(space, 2.0, 1.0);
    }
    
    public P2HardAssociation(Space space, double wellDiameter, double epsilon) {
        super(space);
        setEpsilon(epsilon);
        setWellDiameter(wellDiameter);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
    }
    
   /**
    * Implements the collision dynamics.  Does not deal with the hard cores, only the wells.  This
    * section is essentially the same as PotentialSquareWell without the hard core section.
    */
    public void bump(IAtomList pair, double falseTime) {
        double eps = 1e-6;
        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);

        double reduced_m = 1/(atom0.getType().rm() + atom1.getType().rm());
        double nudge = 0;
        if (bij > 0.0) {    //Separating
            double ke = bij*bij*reduced_m/(2*r2);
            if (ke < epsilon) {    //Not enough energy to escape the well
                lastCollisionVirial = 2*bij*reduced_m;
                nudge = -eps;
                lastEnergyChange = 0.0;
            }
            else {  //Escaping the well
                lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
                nudge = eps;
                lastEnergyChange = epsilon;
            }
        }
        else {          //Approaching
            lastCollisionVirial = reduced_m*(bij + Math.sqrt(bij*bij + 2.0*r2*epsilon/reduced_m));  //might need to double check
            nudge = -eps;
            lastEnergyChange = -epsilon;
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PE(dv);
        atom1.getVelocity().ME(dv);
        atom0.getPosition().Ea1Tv1(-falseTime,dv);
        atom1.getPosition().Ea1Tv1(falseTime,dv);
        if(nudge != 0) {
            atom0.getPosition().PEa1Tv1(-nudge,dr);
            atom1.getPosition().PEa1Tv1(nudge,dr);
        }
    }
    
    
    public double lastCollisionVirial() {return lastCollisionVirial;}
    
    public double energyChange() {return lastEnergyChange;}
    
    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }
    
   /**
    * Computes the next time of collision of the given atomPair assuming free flight.  Only computes the next
    * collision of the wells.  Takes into account both separation and convergence.
    */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double time = Double.POSITIVE_INFINITY;
        
        if (r2 < wellDiameterSquared) {         //check to see if already inside wells
            double discr = bij*bij - v2 * (r2 - wellDiameterSquared);
            time = (-bij + Math.sqrt(discr))/v2;
        }
        else {
            if (bij < 0.) { //Approaching
                double discr = bij*bij - v2 * (r2 - wellDiameterSquared );
                if (discr > 0.) {
                    time = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        return time + falseTime;
    }
    
  /**
   * Returns -epsilon if less than well diameter, or zero otherwise.
   */
    public double u(double r2) {
        return (r2 < wellDiameterSquared) ?  -epsilon : 0.0;
    }
    
   /**
    * Accessor for well diameter.
    * Since the well-diameter is not a multiplier in this potential as in square well, it is necessary
    * to be able to set this manually if so desired.
    */
    public double getWellDiameter() {return wellDiameter;}
    
   /**
    * Accessor for well-diameter.
    * Allows manual changing of the well diameter since it is not merely a multiple of the core-diameter
    * as in square well.
    */
    
    public void setWellDiameter(double w) {
        wellDiameter = w;
        wellDiameterSquared = w*w;
    }
    
    public Dimension getWellDiameterDimension() {
        return Length.DIMENSION;
    }
    
    /**
     * Returns the well diameter.
     */
    public double getRange() {
        return wellDiameter;
    }
    
   /**
    * Accessor method for depth of well.
    */
    public double getEpsilon() {return epsilon;}
    
   /**
    * Accessor method for depth of well.
    */
    public void setEpsilon(double s) {
        epsilon = s;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    private static final long serialVersionUID = 1L;
    private double wellDiameter, wellDiameterSquared;
    private double epsilon;
    private double lastCollisionVirial, lastCollisionVirialr2;
    private final Tensor lastCollisionVirialTensor;
    private double lastEnergyChange;
    private final Vector dv;
}
