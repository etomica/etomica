package etomica;

/**
 * Ideal-gas potential, which defines all atoms to have no interaction
 * Returns a zero energy and force, and an infinite collision time.
 *
 * @author David Kofke
 */
public class P2IdealGas extends Potential2 implements Potential2Hard, Potential2Soft, EtomicaElement {
    
    public String getVersion() {return "P2IdealGas:01.07.07/"+Potential2.VERSION;}

    private final Space.Vector zero;
    private final Space.Tensor zilch;
    
    public P2IdealGas() {
        this(Simulation.instance);
    }
    
    public P2IdealGas(Simulation sim) {
        super(sim);
        zero = sim.space().makeVector();
        zero.E(0.0);
        zilch = sim.space().makeTensor();
        zilch.E(0.0);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Potential corresponding to no intermolecular interactions");
        return info;
    }

   /**
    * Always returns zero
    */
    public double energy(AtomPair pair) {return 0.0;}
    
   /**
    * Always returns false
    */
    public boolean overlap(AtomPair pair) {return false;}

   /**
    * Returns without performing any action to pair
    */
    public void bump(AtomPair pair) {}
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
    public Space.Vector gradient(AtomPair pair) {return zero;}
    
   /**
    * Always returns zero
    */
    public double virial(AtomPair pair) {return 0.0;}
    
   /**
    * Always returns zero
    */
    public double hyperVirial(AtomPair pair) {return 0.0;}
    
    /**
     * Always return zero.
     */
    public double integral(double rc) {return 0.0;}
}

