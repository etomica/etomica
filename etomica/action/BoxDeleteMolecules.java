package etomica.action;

import etomica.atom.AtomFilter;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomIteratorAllMolecules;

/**
 * Deletes molecules from a box as determined by an AtomFilter. Atoms deleted
 * are those for which the filter's accept method returns false.
 * 
 * @author David Kofke
 *  
 */
public class BoxDeleteMolecules extends BoxActionAdapter {

    /**
     * @param filter
     *            determines the atoms that will be deleted by the action; those
     *            for which filter.accept returns false are deleted
     */
    public BoxDeleteMolecules(AtomFilter filter) {
        this.filter = filter;
        iterator = new AtomIteratorAllMolecules();
    }

    /**
     * Performs the action of deleting accept == false molecules, considering
     * all molecules in the box last given to setBox. If no box was given,
     * no action is performed and method returns quietly.
     */
    public void actionPerformed() {
        iterator.setBox(box);
        iterator.reset();
        for (IMolecule molecule = (IMolecule)iterator.nextAtom(); molecule != null;
             molecule = (IMolecule)iterator.nextAtom()) {
            if (!filter.accept(molecule)) {
                box.removeMolecule(molecule);
            }
        }
    }

    private static final long serialVersionUID = 1L;
    private final AtomFilter filter;
    private final AtomIteratorAllMolecules iterator;
}