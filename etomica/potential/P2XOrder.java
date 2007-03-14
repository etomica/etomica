package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Hard potential that enforces ordering of the x-coordinates of the
 * pairs.  Returns infinite energy if the difference in atom indexes
 * and difference in x coordinates are of opposite sign; returns
 * zero otherwise.  Designed for use in 1D simulations.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public class P2XOrder extends Potential2 implements PotentialHard {
    
    private static final long serialVersionUID = 1L;
    protected final IVector dr, drOld;
    protected Phase phase;
    
    public P2XOrder(Space space) {
        super(space);
        dr = space.makeVector();
        drOld = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Potential that enforces ordering in x coordinate (meant for 1D simulations)");
        return info;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(AtomSet pair, double falseTime) {
        throw new RuntimeException("P2XOrder.collisionTime not implemented");
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(AtomSet pair, double falseTime) {
        throw new RuntimeException("P2XOrder.bump not implemented");
    }
    
    public double lastCollisionVirial() {
        throw new RuntimeException("P2XOrder.lastCollisionVirial not implemented");
    }
    
    public Tensor lastCollisionVirialTensor() {
        throw new RuntimeException("P2XOrder.lastCollisionVirialTensor not implemented");
    }
    
    /**
     * Interaction energy of the pair.
     * Zero if x coordinates are ordered differently from atom indexes.
     */
    public double energy(AtomSet pair) {
        dr.Ev1Mv2(((AtomLeaf)((AtomPair)pair).atom1).getCoord().getPosition(), ((AtomLeaf)((AtomPair)pair).atom0).getCoord().getPosition());
        drOld.E(dr);
        int dI = ((AtomPair)pair).atom1.getIndex() - ((AtomPair)pair).atom0.getIndex();
        if (dI == ((AtomPair)pair).atom1.parentGroup().childAtomCount()-1) {
            dr.ME(phase.getBoundary().getDimensions());
            return (dr.x(0) * dI > 0.0) ? Double.POSITIVE_INFINITY : 0.0;
        }
        else if (dI == 1) {
            return (dr.x(0) * dI < 0.0) ? Double.POSITIVE_INFINITY : 0.0;
        }
        else {
            return 0;
        }
    }
    
    public double energyChange() {
        throw new RuntimeException("P2XOrder.energyChange not implemented");
    }
    
    /**
     * Returns infinity.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;
    }
    
}//end of P2XOrder