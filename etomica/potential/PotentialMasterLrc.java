package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.space.Space;

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

    protected PotentialMasterLrc(Space space) {
        super(space);
    }
    
    /**
     * Performs given PotentialCalculation on all LRC potentials added to this group.
     * Checks that group is enabled, box is not null, that it has lrcEnabled,
     * and that the given IteratorDirective has includeLrc set to true; if all
     * are so, calculation is performed.
     */
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled || box == null || !box.isLrcEnabled() || !id.includeLrc) return;
        IAtom targetAtom = id.getTargetAtom();
        boolean boxChanged = (box != mostRecentBox);
        mostRecentBox = box;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(!link.enabled) continue;
            if(boxChanged) {
                link.iterator.setBox(box);
                link.potential.setBox(box);
            }
            link.iterator.setTarget(targetAtom);
            ((Potential0Lrc)link.potential).setTargetAtoms(targetAtom);
            pc.doCalculation(link.iterator, id, link.potential);
        }
    }
    
    private static final long serialVersionUID = 1L;
}
