package etomica;

/**
 * Hard potential that enforces ordering of the x-coordinates of the
 * pairs.  Returns infinite energy if the difference in atom indexes
 * and difference in x coordinates are of opposite sign; returns
 * zero otherwise.  Designed for use in 1D simulations.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public class P2XOrder extends Potential2Hard implements EtomicaElement {
    
   public String getVersion() {return "P2XOrder:02.04.15/"+Potential2.VERSION;}
   
   protected final Space.Vector dr;
    
    public P2XOrder() {
        this(Simulation.instance.hamiltonian.potential);
    }

    public P2XOrder(PotentialGroup parent) {
        super(parent);
        dr = parentSimulation().space().makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Potential that enforces ordering in x coordinate (meant for 1D simulations)");
        return info;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(AtomPair pair) {
        throw new RuntimeException("P2XOrder.collisionTime not implemented");
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(AtomPair pair) {
        throw new RuntimeException("P2XOrder.bump not implemented");
    }
    
    public double lastCollisionVirial() {
        throw new RuntimeException("P2XOrder.lastCollisionVirial not implemented");
    }
    
    public Space.Tensor lastCollisionVirialTensor() {
        throw new RuntimeException("P2XOrder.lastCollisionVirialTensor not implemented");
    }
    
    /**
     * Interaction energy of the pair.
     * Zero if x coordinates are ordered differently from atom indexes.
     */
    public double energy(AtomPair pair) {
 //       double deltaX = pair.dr(0);
        double deltaX = pair.atom2().coord.position(0) - pair.atom1().coord.position(0);
        int dI = pair.atom2().node.index() - pair.atom1().node.index();
        return (deltaX * dI < 0.0) ? Double.MAX_VALUE : 0.0;
    }
    
    
}//end of P2XOrder