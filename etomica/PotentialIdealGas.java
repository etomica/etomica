package simulate;

/**
 * Ideal-gas potential, which defines all atoms to have no interaction
 * Returns a zero energy and force, and an infinite collision time
 */
public class PotentialIdealGas extends Potential implements Potential.Hard, Potential.Soft {
    
    private final Space.Vector zero;
    private final Space.Tensor zilch;
    
    public PotentialIdealGas() {
        this(Simulation.instance);
    }
    
    public PotentialIdealGas(Simulation sim) {
        super(sim);
        zero = sim.space().makeVector();
        zero.E(0.0);
        zilch = sim.space().makeTensor();
        zilch.E(0.0);
    }
    
   /**
    * Always returns zero
    */
    public double energy(AtomPair pair) {return 0.0;}
   /**
    * Always returns zero
    */
    public double energyLRC(int n1, int n2, double V) {return 0.0;}
   /**
    * Always returns zero
    */
    public double pressureLRC(int n1, int n2, double V) {return 0.0;}
    
   /**
    * Always returns false
    */
    public boolean overlap(AtomPair pair) {return false;}

   /**
    * Returns without performing any action to pair
    */
    public void bump(AtomPair pair) {return;}
   /**
    * Always returns infinity (Double.MAX_VALUE)
    */
    public double collisionTime(AtomPair pair) {return Double.MAX_VALUE;}
   /**
    * Always returns zero
    */
    public double lastCollisionVirial() {return 0.0;}
   /**
    * Always returns a zero tensor
    */
    public Space.Tensor lastCollisionVirialTensor() {return zilch;}
   /**
    * Always returns a zero vector
    */
    public Space.Vector force(AtomPair pair) {return zero;}
    
   /**
    * Always returns zero
    */
    public double virial(AtomPair pair) {return 0.0;}
    
}

