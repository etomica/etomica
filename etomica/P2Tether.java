package etomica;
import etomica.units.Dimension;
/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction.
 * Interaction of atoms is zero if separated by less than the tether length.  Atoms
 * undergo an impulsive attractive collision when attempting to separate by more than the tether distance.
 *
 * @author David Kofke
 */
public class P2Tether extends Potential2Hard implements EtomicaElement {

  public String getVersion() {return "P2HardSphere:01.07.03/"+Potential2.VERSION;}

  private double tetherLength, tetherLengthSquared;
  private double lastCollisionVirial = 0.0;
  private double lastCollisionVirialr2 = 0.0;
  private final Space.Vector dr;
  private final Space.Tensor lastCollisionVirialTensor;

  public P2Tether() {
    this(Simulation.instance);
  }
  public P2Tether(Simulation sim) {
    super(sim);
    setTetherLength(0.75*Default.ATOM_SIZE);
    lastCollisionVirialTensor = sim.space().makeTensor();
    dr = sim.space().makeVector();
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
  public final void bump(AtomPair pair) {
        double r2 = pair.r2();
        dr.E(pair.dr());
        lastCollisionVirial = 2.0/(pair.atom1().coord.rm() + pair.atom2().coord.rm())*pair.vDotr();
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        pair.cPair.push(lastCollisionVirialr2);
  }


    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    public final Space.Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
    

  
  /**
   * Time at which two atoms will reach the end of their tether, assuming free-flight kinematics
   */
  public final double collisionTime(AtomPair pair) {
    boolean minus;  //flag for +/- root
    double r2 = pair.r2();
    double bij = pair.vDotr();
 //       if(r2 < sig2) {return (bij > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
    if(r2 > tetherLengthSquared && bij >= 0) {return 0.0;}  //outside tether, moving apart; collide now
    double v2 = pair.v2();
    double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
    return (-bij + Math.sqrt(discr))/v2;
    }
  
  /**
   * Returns infinity if separation is greater than tether distance, zero otherwise
   */
    public double energy(AtomPair pair) {
        return (pair.r2() > tetherLengthSquared) ? Double.MAX_VALUE : 0.0;
    }
      
}//end of P2Tether
  