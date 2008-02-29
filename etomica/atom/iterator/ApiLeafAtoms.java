package etomica.atom.iterator;

import etomica.api.IBox;

/**
 * Iterator that returns all pairs that can be formed from all leaf atoms of a
 * given box. Wraps an ApiIntraList instance.
 */
public class ApiLeafAtoms extends AtomsetIteratorAdapter implements
        AtomsetIteratorBoxDependent {

    /**
     * Creates new pair iterator that requires reset() before beginning
     * iteration.
     */
    public ApiLeafAtoms() {
        super(new ApiIntraArrayList());
    }

    /**
     * Conditions iterator to return all leaf-atom pairs from the given box.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(IBox box) {
        ((ApiIntraArrayList)iterator).setList(box.getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
