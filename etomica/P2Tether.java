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
public class P2Tether extends Potential2 implements PotentialHard {

  private double tetherLength, tetherLengthSquared;
  private double lastCollisionVirial = 0.0;
  private double lastCollisionVirialr2 = 0.0;
  private final Space.Vector dr;
  private final Space.Tensor lastCollisionVirialTensor;

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
  public final void bump(Atom[] pair, double falseTime) {
  		cPair.trueReset(pair[0].coord,pair[1].coord,falseTime);
        double r2 = cPair.r2();
        dr.E(cPair.dr());
        lastCollisionVirial = 2.0/(pair[0].coord.rm() + pair[1].coord.rm())*cPair.vDotr();
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dr.TE(lastCollisionVirialr2);
        cPair.truePush(dr,falseTime);
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
  public final double collisionTime(Atom[] pair, double falseTime) {
  	cPair.trueReset(pair[0].coord,pair[1].coord,falseTime);
    double r2 = cPair.r2();
    double bij = cPair.vDotr();
    if(Default.FIX_OVERLAP && r2 > tetherLengthSquared && bij > 0) {return 0.0;}  //outside tether, moving apart; collide now
    double v2 = cPair.v2();
    double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
    return (-bij + Math.sqrt(discr))/v2 + falseTime;
  }
  
  /**
   * Returns infinity if separation is greater than tether distance, zero otherwise
   */
    public double energy(Atom[] pair) {
//        System.out.println(pair.atom1.toString()+"  "+pair.atom2.toString()+"   "+pair.r2()+"   "+tetherLengthSquared);
    	cPair.reset(pair[0].coord,pair[1].coord);
        return (cPair.r2() > tetherLengthSquared) ? Double.MAX_VALUE : 0.0;
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
  