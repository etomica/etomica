package etomica.potential;
import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.space.CoordinatePairKinetic;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Dimension;
/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction.
 * Interaction of atoms is zero if separated by less than the tether length.  Atoms
 * undergo an impulsive attractive collision when attempting to separate by more than the tether distance.
 *
 * @author David Kofke
 */
public class P2Tether extends Potential2HardSpherical {

  private double tetherLength, tetherLengthSquared;
  private double lastCollisionVirial = 0.0;
  private double lastCollisionVirialr2 = 0.0;
  private final Vector dr;
  private final Tensor lastCollisionVirialTensor;

  public P2Tether() {
    this(Simulation.getDefault().space);
  }
  public P2Tether(Space space) {
    super(space);
    setTetherLength(0.75*Default.ATOM_SIZE);
    lastCollisionVirialTensor = space.makeTensor();
    dr = space.makeVector();
  }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard string between adjacent atoms, hard sphere for non-adjacents");
        return info;
    }

  /**
   * Accessor method for the tether distance
   */
  public final double getTetherLength() {return tetherLength;}
  /**
   * Accessor method for the tether distance
   */
  public final void setTetherLength(double t) {
      tetherLength = t;
      tetherLengthSquared = t*t;
  }
  public final Dimension getTetherLengthDimension() {return Dimension.LENGTH;}

  /**
   * Implements collision dynamics for pair attempting to separate beyond tether distance
   */
  public final void bump(AtomSet pair, double falseTime) {
      cPair.reset((AtomPair)pair);
      ((CoordinatePairKinetic)cPair).resetV();
      dr.E(cPair.dr());
      Vector dv = ((CoordinatePairKinetic)cPair).dv();
      dr.PEa1Tv1(falseTime,dv);
      double r2 = dr.squared();
      double bij = dr.dot(dv);
        lastCollisionVirial = 2.0/(((AtomPair)pair).atom0.type.rm() + ((AtomPair)pair).atom1.type.rm())*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        ((ICoordinateKinetic)((AtomPair)pair).atom0.coord).velocity().PE(dv);
        ((ICoordinateKinetic)((AtomPair)pair).atom1.coord).velocity().ME(dv);
        ((AtomPair)pair).atom0.coord.position().Ea1Tv1(-falseTime,dv);
        ((AtomPair)pair).atom1.coord.position().Ea1Tv1(falseTime,dv);
  }


    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
    

  
  /**
   * Time at which two atoms will reach the end of their tether, assuming free-flight kinematics
   */
  public final double collisionTime(AtomSet pair, double falseTime) {
      cPairNbr.reset((AtomPair)pair);
      ((CoordinatePairKinetic)cPairNbr).resetV();
      dr.E(cPairNbr.dr());
      Vector dv = ((CoordinatePairKinetic)cPairNbr).dv();
      dr.Ea1Tv1(falseTime,dv);
      double r2 = dr.squared();
      double bij = dr.dot(dv);
      double v2 = dv.squared();
      if(Default.FIX_OVERLAP && r2 > tetherLengthSquared && bij > 0) {return 0.0;}  //outside tether, moving apart; collide now
      double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
      return (-bij + Math.sqrt(discr))/v2 + falseTime;
  }
  
  /**
   * Returns infinity if separation is greater than tether distance, zero otherwise
   */
    public double u(double r2) {
        return (r2 > tetherLengthSquared) ? Double.MAX_VALUE : 0.0;
    }
    
    public double energyChange() {return 0.0;}
    
    /**
     * Method to demonstrate use of this potential is given as the
     * main method of the SpeciesSpheres class.
     */
/*    public static void main(String[] args) {
        SpeciesSpheres.main(args);
    }
 */     
}//end of P2Tether
  