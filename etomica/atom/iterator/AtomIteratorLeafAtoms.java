package etomica.atom.iterator;

import etomica.box.Box;

/**
 * Iterator that will loop over all leaf atoms in a box. Can be configured to
 * iterate all leaf atoms, or only those of a particular species.
 */
public final class AtomIteratorLeafAtoms extends AtomIteratorAdapter implements 
        AtomIteratorBoxDependent {

    /**
     * Creates iterator with no box specified. Iteration will return no atoms
     * until a call to setBox is performed.
     */
    public AtomIteratorLeafAtoms() {
        super(new AtomIteratorArrayListSimple());
    }

    /**
     * Creates iterator conditioned to give all leaf atoms of the specified
     * box. Call to reset() is required before beginning iteration.
     */
    public AtomIteratorLeafAtoms(Box box) {
        this();
        setBox(box);
    }

    /**
     * Configures iterator to form its iterates from the leaf atoms of the given
     * box. If a species was previously (or subsequently) set, iterates will
     * be the leaf atoms of under the species in the specified box.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        ((AtomIteratorArrayListSimple)iterator).setList(box.getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
