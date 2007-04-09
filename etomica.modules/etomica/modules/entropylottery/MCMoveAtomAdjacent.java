package etomica.modules.entropylottery;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.phase.Phase;
import etomica.space.BoundaryPeriodic;
import etomica.space.IVector;
import etomica.util.IRandom;

/**
 * Monte Carlo move that moves an atom by +/- 1 unit in a random dimension.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomAdjacent extends MCMovePhase {
    
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected IVector translationVector;
    protected Atom atom;
    protected AtomSource atomSource;
    protected final IRandom random;

    public MCMoveAtomAdjacent(IRandom random) {
        super(null);
        this.random = random;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        perParticleFrequency = true;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        translationVector.E(0);
        int i = random.nextInt(translationVector.getD());
        translationVector.setX(i, random.nextInt(2)*2-1);
        ((AtomLeaf)atom).getCoord().getPosition().PE(translationVector);
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
        IVector position = ((AtomLeaf)atom).getCoord().getPosition();
        IVector dimensions = phase.getBoundary().getDimensions();
        for (int i=0; i<position.getD(); i++) {
            // if we're non-periodic, ensure we didn't try to jump over the boundary
            int x = (int)Math.round(position.x(i)+dimensions.x(i)*0.5-0.5);
            if (x < 0 || x >= (int)Math.round(dimensions.x(i))) {
                if (!periodicity[i]) {
                    // failure
                    return Double.NEGATIVE_INFINITY;
                }
                //wrap around -- OK, it's a hack.  deal with it.
                if (x < 0) {
                    position.setX(i, position.x(i)+dimensions.x(i));
                }
                else {
                    position.setX(i, 0);
                }
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
        ((AtomLeaf)atom).getCoord().getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        translationVector = phase.getSpace().makeVector();
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