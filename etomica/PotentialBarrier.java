package etomica;

/**
 */
public class PotentialBarrier extends Potential implements Potential.Hard {

  protected double barrierDiameter, barrierDiameterSquared;
  protected double barrier;
  protected double lastCollisionVirial, lastCollisionVirialr2;
  protected Space.Tensor lastCollisionVirialTensor;
  protected Space.Vector dr;
   
  public PotentialBarrier() {
    this(Simulation.instance);
  }
  public PotentialBarrier(Simulation sim) {
    this(sim, Default.ATOM_SIZE, Default.POTENTIAL_WELL);
  }
  public PotentialBarrier(double barrierDiameter, double barrier) {
    this(Simulation.instance, barrierDiameter, barrier);
  }
  public PotentialBarrier(Simulation sim, double barrierDiameter, double barrier) {
    super(sim);
    setBarrierDiameter(barrierDiameter);
    setBarrierEnergy(barrier);
    dr = sim.space().makeVector();
    lastCollisionVirialTensor = sim.space().makeTensor();
  }

/**
 * Returns true if separation of pair is less than core diameter, false otherwise
 */
  public boolean overlap(AtomPair pair) {return false;}

/**
 * Implements collision dynamics between two square-barrier atoms.
 * Includes all possibilities involving collision of hard cores, and collision of barriers
 * both approaching and diverging
 */
  public void bump(AtomPair pair) {
    double eps = 1.0e-6;
    double r2 = pair.r2();
    double bij = pair.vDotr();
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/(pair.atom1().rm() + pair.atom2().rm());
    dr.E(pair.dr());
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;

      if(r2 < barrierDiameterSquared) return;
      if(bij > 0.0) {         // Separating
	    if(ke < barrier) {     // Not enough kinetic energy to escape
	       lastCollisionVirial = 2.0*reduced_m*bij;
	       r2New = (1-eps)*barrierDiameterSquared; 
	    }
	    else {                 // Escape
//	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - barrier)))/r;
	       lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*barrier/reduced_m));
	       r2New = (1+eps)*barrierDiameterSquared;
	       if(pair.atom1().atomLink != null) {
	            pair.atom1().atomLink[0] = null;
	            pair.atom2().atomLink[0] = null;
	       }
	    }
      }
      else {// Approaching
        if(ke < barrier) { //not enough kinetic energy to form bond
	       lastCollisionVirial = 2.0*reduced_m*bij;
	       r2New = (1+eps)*barrierDiameterSquared; 
        }
        else {  //bonding
	     lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*barrier/reduced_m));
	     r2New = (1-eps)*barrierDiameterSquared;
	     if(pair.atom1().atomLink != null) {
	        //remove any existing bonds pair atoms have with other atoms
	        if(pair.atom1().atomLink[0] != null) {
	            pair.atom1().atomLink[0].atom().atomLink[0] = null;
	        }
	        if(pair.atom2().atomLink[0] != null) {
	            pair.atom2().atomLink[0].atom().atomLink[0] = null;
	        }
	        //bond pair atoms
	        pair.atom1().atomLink[0] = new Atom.Linker(pair.atom2());
	        pair.atom2().atomLink[0] = new Atom.Linker(pair.atom1());
	     }
	    }
      }

    lastCollisionVirialr2 = lastCollisionVirial/r2;
    pair.cPair.push(lastCollisionVirialr2);
    if(r2New != r2) pair.cPair.setSeparation(r2New);
  }//end of bump method
  
 /**
  * Always returns zero (not yet implemented)
  */
  public double lastCollisionVirial() {
    return lastCollisionVirial;
  }

  public Space.Tensor lastCollisionVirialTensor() {
    lastCollisionVirialTensor.E(dr, dr);
    lastCollisionVirialTensor.TE(lastCollisionVirialr2);
    return lastCollisionVirialTensor;
  }

/**
 * Computes next time of collision of two square-barrier atoms, assuming free-flight kinematics.
 * Collision may occur when cores collides, or when barriers first encounter each other on
 * approach, or when they edge of the barriers are reached as atoms diverge.
 */
  public double collisionTime(AtomPair pair) {
    double discr = 0.0;

    double bij = pair.vDotr();
    double r2 = pair.r2();
    double v2 = pair.v2();

    double tij = Double.MAX_VALUE;

    if(r2 < barrierDiameterSquared) {  // Already inside barriers
        return Double.MAX_VALUE;
//	      discr = bij*bij - v2 * ( r2 - barrierDiameterSquared );  // This is always > 0
//	      tij = (-bij + Math.sqrt(discr))/v2;
    }
    else {              // Outside barriers; look for collision at barrier
        if(bij < 0.0) {
	        discr = bij*bij - v2 * ( r2 - barrierDiameterSquared );
	        if(discr > 0) {
    	        tij = (-bij - Math.sqrt(discr))/v2;
	        }
        }
    }
    return tij;
  }
  
  /**
   * Returns infinity if overlapping, -barrier if otherwise less than barrier diameter, or zero if neither.
   */
    public double energy(AtomPair pair) {
        return 0.0;
    }
    
    /**
     * Always returns zero.
     */
    public double energyLRC(int n1, int n2, double V) {return 0.0;}
 
    /**
     * Accessor method for barrier diameter.
     */
    public double getBarrierDiameter() {return barrierDiameter;}
    /**
     * Accessor method for barrier diameter.
     */
    public void setBarrierDiameter(double c) {
        barrierDiameter = c;
        barrierDiameterSquared = c*c;
    }

   /**
    * Accessor method for depth of barrier
    */
    public double getBarrierEnergy() {return barrier;}
   /**
    * Accessor method for depth of barrier
    */
    public void setBarrierEnergy(double eps) {
        barrier = eps;
    }

}
  