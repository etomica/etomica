package etomica.modules.entropylottery;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryPeriodic;
import etomica.space.Vector;

/**
 * Monte Carlo move that moves an atom by +/- 1 unit in a random dimension.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomAdjacent extends MCMovePhase {
    
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected Vector translationVector;
    protected Atom atom;
    protected AtomSource atomSource;

    public MCMoveAtomAdjacent() {
        super(null);
        atomSource = new AtomSourceRandomLeaf();
        perParticleFrequency = true;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        translationVector.E(0);
        int i = Simulation.random.nextInt(translationVector.D());
        translationVector.setX(i, Simulation.random.nextInt(2)*2-1);
        ((AtomLeaf)atom).coord.position().PE(translationVector);
        return true;
    }
    
    public double getA() {return 1.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        boolean[] periodicity = ((BoundaryPeriodic)phase.getBoundary()).getPeriodicity();
        Vector position = ((AtomLeaf)atom).coord.position();
        Vector dimensions = phase.getBoundary().getDimensions();
        for (int i=0; i<position.D(); i++) {
            if (periodicity[i]) continue;
            // if we're non-periodic, ensure we didn't try to jump over the boundary
            int x = (int)Math.round(position.x(i)+phase.getBoundary().getDimensions().x(i)*0.5);
            if (x < 0 || x > (int)Math.round(dimensions.x(i)-1)) {
                // failure
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0;
    }
    
    public double energyChange() {return 0;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        translationVector.TE(-1);
        ((AtomLeaf)atom).coord.position().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        translationVector = phase.space().makeVector();
        atomSource.setPhase(p);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }

    private static final long serialVersionUID = 1L;
}