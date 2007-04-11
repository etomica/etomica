package etomica.action;

import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.iterator.AtomIteratorAllMolecules;

/**
 * Deletes molecules from a phase as determined by an AtomFilter. Atoms deleted
 * are those for which the filter's accept method returns false.
 * 
 * @author David Kofke
 *  
 */
public class PhaseDeleteMolecules extends PhaseActionAdapter {

    /**
     * @param filter
     *            determines the atoms that will be deleted by the action; those
     *            for which filter.accept returns false are deleted
     */
    public PhaseDeleteMolecules(AtomFilter filter) {
        this.filter = filter;
        iterator = new AtomIteratorAllMolecules();
    }

    /**
     * Performs the action of deleting accept == false molecules, considering
     * all molecules in the phase last given to setPhase. If no phase was given,
     * no action is performed and method returns quietly.
     */
    public void actionPerformed() {
        iterator.setPhase(phase);
        iterator.reset();
        while (iterator.hasNext()) {
            Atom molecule = iterator.nextAtom();
            if (!filter.accept(molecule)) {
                molecule.getParentGroup().removeChildAtom(molecule);
            }
        }
    }

    private static final long serialVersionUID = 1L;
    private final AtomFilter filter;
    private final AtomIteratorAllMolecules iterator;
}