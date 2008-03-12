package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.IteratorDirective;
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
    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled || !id.includeLrc) return;
        IAtom targetAtom = id.getTargetAtom();
        boolean boxChanged = (box != mostRecentBox);
        mostRecentBox = box;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(!link.enabled) continue;
            final IPotential potential = link.potential;
            final AtomsetIteratorPDT atomIterator = link.iterator;
            if(boxChanged) {
                atomIterator.setBox(box);
                potential.setBox(box);
            }
            atomIterator.setTarget(targetAtom);
            ((Potential0Lrc)potential).setTargetAtoms(targetAtom);
            if (potential instanceof PotentialGroup) {
                ((PotentialGroup)potential).calculate(atomIterator, id, pc);
            }
            else {
                atomIterator.reset();
                for (IAtomSet atoms = atomIterator.next(); atoms != null;
                     atoms = atomIterator.next()) {
                    pc.doCalculation(atoms, potential);
                }
            }
        }
    }
    
    private static final long serialVersionUID = 1L;
}
