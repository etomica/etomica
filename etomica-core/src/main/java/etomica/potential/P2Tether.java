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
import etomica.units.Length;
/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction.
 * Interaction of atoms is zero if separated by less than the tether length.  Atoms
 * undergo an impulsive attractive collision when attempting to separate by more than the tether distance.
 *
 * @author David Kofke
 */
public class P2Tether extends Potential2HardSpherical {

    public P2Tether(Space space) {
        this(space, 0.75, false);
    }
    
    public P2Tether(Space space, double tetherLength, boolean ignoreOverlap) {
        super(space);
        setTetherLength(tetherLength);
        this.ignoreOverlap = ignoreOverlap;
        lastCollisionVirialTensor = space.makeTensor();
        dv = space.makeVector();
    }

    /**
     * Accessor method for the tether distance
     */
    public double getTetherLength() {return tetherLength;}
    /**
     * Accessor method for the tether distance
     */
    public void setTetherLength(double t) {
        tetherLength = t;
        tetherLengthSquared = t*t;
    }
    public Dimension getTetherLengthDimension() {return Length.DIMENSION;}

    /**
     * Implements collision dynamics for pair attempting to separate beyond tether distance
     */
    public final void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        lastCollisionVirial = 2.0/(rm0 + rm1)*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1( rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);
    }


    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }

    /**
     * Returns the tether length. The potential range is in fact infinite, but if
     * the integrator is generating configurations correctly, there will be no atoms
     * interacting beyond the tetherlength distance.  
     */
    public double getRange() {
        return tetherLength;
    }
  
    /**
     * Time at which two atoms will reach the end of their tether, assuming free-flight kinematics
     */
    public final double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic coord0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic coord1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        if(ignoreOverlap && r2 > tetherLengthSquared && bij > 0) {return 0.0;}  //outside tether, moving apart; collide now
        double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
        return (-bij + Math.sqrt(discr))/v2 + falseTime;
    }
  
    /**
     * Returns infinity if separation is greater than tether distance, zero otherwise
     */
    public double u(double r2) {
        return (r2 > tetherLengthSquared) ? Double.POSITIVE_INFINITY : 0.0;
    }
    
    public double energyChange() {return 0.0;}
    
    private static final long serialVersionUID = 1L;
    private double tetherLength, tetherLengthSquared;
    private double lastCollisionVirial = 0.0;
    private double lastCollisionVirialr2 = 0.0;
    private final Vector dv;
    private final Tensor lastCollisionVirialTensor;
    private final boolean ignoreOverlap;
   
}//end of P2Tether
