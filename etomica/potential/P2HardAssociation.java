package etomica.potential;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Purely attractive square-well potential with no repulsive core.  Similar
 * to P2Tether, but permits atoms to separate.
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
public class P2HardAssociation extends Potential2HardSpherical {

    public P2HardAssociation(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().potentialCutoffFactor, sim.getDefaults().potentialWell);
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
    public void bump(AtomSet pair, double falseTime) {
        double eps = 1e-6;
        AtomLeaf atom0 = (AtomLeaf)((AtomPair)pair).atom0;
        AtomLeaf atom1 = (AtomLeaf)((AtomPair)pair).atom1;
        ICoordinateKinetic coord0 = (ICoordinateKinetic)atom0.getCoord();
        ICoordinateKinetic coord1 = (ICoordinateKinetic)atom1.getCoord();
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        nearestImageTransformer.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);

        double reduced_m = 1/(((AtomTypeLeaf)atom0.getType()).rm() + ((AtomTypeLeaf)atom1.getType()).rm());
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
        coord0.getVelocity().PE(dv);
        coord1.getVelocity().ME(dv);
        coord0.getPosition().Ea1Tv1(-falseTime,dv);
        coord1.getPosition().Ea1Tv1(falseTime,dv);
        if(nudge != 0) {
            coord0.getPosition().PEa1Tv1(-nudge,dr);
            coord1.getPosition().PEa1Tv1(nudge,dr);
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
    public double collisionTime(AtomSet pair, double falseTime) {
        ICoordinateKinetic coord0 = (ICoordinateKinetic)((AtomLeaf)((AtomPair)pair).atom0).getCoord();
        ICoordinateKinetic coord1 = (ICoordinateKinetic)((AtomLeaf)((AtomPair)pair).atom1).getCoord();
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        nearestImageTransformer.nearestImage(dr);

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
    private final IVector dv;
}
