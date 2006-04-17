package etomica.integrator.mcmove;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;

/**
 * MCMove parent class for moves that operate on a single Phase.
 *
 * @author Andrew Schultz
 */
public abstract class MCMovePhase extends MCMove {

    public MCMovePhase(PotentialMaster potentialMaster) {
        this(potentialMaster, new MCMoveTracker());
    }

    public MCMovePhase(PotentialMaster potentialMaster,
            MCMoveTracker acceptanceTracker) {
        super(potentialMaster, acceptanceTracker);
        perParticleFrequency = false;
    }

    /**
     * Sets the phase on which this move acts.  The phase itself can be changed
     * via this method, if desired.
     */
    public void setPhase(Phase p) {
        phase = p;
    }

    /**
     * @return the phase(s) on which this move acts.
     */
    public Phase getPhase() {
        return phase;
    }

    public AtomIterator affectedAtoms(Phase aPhase) {
        if (phase == aPhase) {
            return affectedAtoms();
        }
        return AtomIteratorNull.INSTANCE;
    }

    public abstract AtomIterator affectedAtoms();
    
    public double energyChange(Phase aPhase) {
        if (phase == aPhase) {
            return energyChange();
        }
        return 0;
    }
    
    public abstract double energyChange();

    /**
     * Indicates whether this move should nominally be performed at a frequency
     * proportional to the number of molecules in the phase.
     * 
     * @see #nominalFrequency
     */
    public final boolean isNominallyPerParticleFrequency() {
        return perParticleFrequency;
    }

    /**
     * Flag indicating whether nominal frequency is interpreted as a
     * perParticleFrequency, or as a full frequency. Default is false, but may
     * be given a different value by subclasses.
     */
    protected boolean perParticleFrequency;
    
    protected Phase phase;
}
