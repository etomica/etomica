package etomica;

/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction.
 * Interaction of atoms is zero if separated by less than the tether length.  Atoms
 * undergo an impulsive attractive collision when attempting to separate by more than the tether distance.
 */
public class PotentialTether extends Potential implements Potential.Hard, EtomicaElement {

  private double tetherLength, tetherLengthSquared;
  private double lastCollisionVirial = 0.0;
  private double lastCollisionVirialr2 = 0.0;
  private final Space.Vector dr;
  private final Space.Tensor lastCollisionVirialTensor;

  public PotentialTether() {
    this(Simulation.instance);
  }
  public PotentialTether(Simulation sim) {
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
   * Always returns false
   */
  public boolean overlap(AtomPair pair) {return false;}

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

  /**
   * Implements collision dynamics for pair attempting to separate beyond tether distance
   */
  public final void bump(AtomPair pair) {
        double r2 = pair.r2();
        dr.E(pair.dr());
        lastCollisionVirial = 2.0/(pair.atom1().rm() + pair.atom2().rm())*pair.vDotr();
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
  
  /**
   * Always returns zero
   */
    public double energyLRC(int n1, int n2, double V) {return 0.0;}
      
}
  