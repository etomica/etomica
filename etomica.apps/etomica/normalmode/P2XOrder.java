package etomica.normalmode;

import etomica.EtomicaInfo;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.Potential2Spherical;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Hard potential that enforces ordering of the x-coordinates of the
 * pairs.  Returns infinite energy if the difference in atom indexes
 * and difference in x coordinates are of opposite sign; returns
 * zero otherwise.  Designed for use in 1D simulations.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public class P2XOrder extends Potential2 {
    
    private static final long serialVersionUID = 1L;
    protected final IVector dr;
    protected Phase phase;
    protected Potential2Spherical wrappedPotential;
    
    public P2XOrder(Space space, Potential2Spherical wrappedPotential) {
        super(space);
        dr = space.makeVector();
        this.wrappedPotential = wrappedPotential;
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Potential that enforces ordering in x coordinate (meant for 1D simulations)");
        return info;
    }

    /**
     * Interaction energy of the pair.
     * Zero if x coordinates are ordered differently from atom indexes.
     */
    public double energy(AtomSet pair) {
        dr.Ev1Mv2(((IAtomPositioned)((AtomPair)pair).atom1).getPosition(), ((IAtomPositioned)((AtomPair)pair).atom0).getPosition());
        int dI = ((AtomPair)pair).atom1.getIndex() - ((AtomPair)pair).atom0.getIndex();
        if (Math.abs(dI) == ((AtomPair)pair).atom1.getParentGroup().getChildList().size()-1) {
            phase.getBoundary().nearestImage(dr);
            return (dr.x(0) * dI > 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else if (dI == 1 || dI == -1) {
            return (dr.x(0) * dI < 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else {
            return 0;
        }
    }
    
    /**
     * Returns infinity.
     */
    public double getRange() {
        return ((Potential)wrappedPotential).getRange();
    }
    
    public Potential2Spherical getWrappedPotential() {
        return wrappedPotential;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;
        ((Potential)wrappedPotential).setPhase(newPhase);
    }
    
}//end of P2XOrder