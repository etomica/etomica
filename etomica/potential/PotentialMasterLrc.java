package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.iterator.IteratorDirective;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

/**
 * Collects potentials used for long-range correction.
 * One instance of this class is added to the PotentialMaster of
 * a Simulation when the first Potential0Lrc class is instantiated.
 * All subsequently created Potential0Lrc classes are added to this instance.
 *
 * @see Potential0Lrc
 * @author David Kofke
 */
 
public class PotentialMasterLrc extends PotentialMaster {

    protected PotentialMasterLrc(Simulation sim) {
        super(sim);
    }
    
    /**
     * Performs given PotentialCalculation on all LRC potentials added to this group.
     * Checks that group is enabled, phase is not null, that it has lrcEnabled,
     * and that the given IteratorDirective has includeLrc set to true; if all
     * are so, calculation is performed.
     */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled || phase == null || !phase.isLrcEnabled() || !id.includeLrc) return;
        IAtom targetAtom = id.getTargetAtom();
        boolean phaseChanged = (phase != mostRecentPhase);
        mostRecentPhase = phase;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(!link.enabled) continue;
            if(phaseChanged) {
                link.iterator.setPhase(phase);
                link.potential.setPhase(phase);
            }
            link.iterator.setTarget(targetAtom);
            ((Potential0Lrc)link.potential).setTargetAtoms(targetAtom);
            pc.doCalculation(link.iterator, id, link.potential);
        }
    }
    
    private static final long serialVersionUID = 1L;
}
