package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Returns all pairs formed from two different untabbed lists of atoms.
 * Incorrect behavior will result if both lists refer to the same instance.
 *  
 */
public class ApiInterArrayList implements AtomsetIterator, java.io.Serializable {

    /**
     * Construct iterator with an empty lists. No iterates will be given until
     * non-empty lists are specified via setList.
     */
    public ApiInterArrayList() {
        this(new AtomArrayList(), new AtomArrayList());
    }

    /**
     * Constructs iterator to return pairs from the given lists. Requires
     * reset() before first use.
     * 
     * @throws IllegalArgumentException
     *             if both lists refer to the same instance
     */
    public ApiInterArrayList(AtomSet outerList, AtomSet innerList) {
        if (outerList == innerList) {
            throw new IllegalArgumentException(
                    "ApiInterList will not work if both iterators are the same instance");
        }
        setOuterList(outerList);
        setInnerList(innerList);
    }

    /**
     * Sets iterator in condition to begin iteration.
     * 
     * @throws IllegalStateException
     *             if outer and inner lists have been set to the same instance
     */
    public void reset() {
        if (outerList == innerList) {
            throw new IllegalStateException(
                    "ApiInterList will not work correctly if inner and outer lists are the same instance");
        }
        outerIndex = 0;
        if (outerList.getAtomCount() == 0) {
            innerIndex = innerList.getAtomCount() - 1;
            return;
        }
        innerIndex = -1;
        atoms.atom0 = outerList.getAtom(outerIndex);
    }

    /**
     * Sets iterator such that hasNext is false.
     */
    public void unset() {
        if (outerList != null) {
            outerIndex = outerList.getAtomCount() - 1;
        }
        if (innerList != null) {
            innerIndex = innerList.getAtomCount() - 1;
        }
    }

    /**
     * Returns the next iterate pair. Returns null if there are no more
     * iterates.
     */
    public AtomSet next() {
        if (innerIndex > innerList.getAtomCount() - 2) {
            if (outerIndex > outerList.getAtomCount() - 2 || innerList.getAtomCount() == 0) {
                return null;
            }
            outerIndex++;
            atoms.atom0 = outerList.getAtom(outerIndex);
            innerIndex = -1;
        }
        innerIndex++;
        atoms.atom1 = innerList.getAtom(innerIndex);
        return atoms;
    }

    /**
     * Performs given action on all pairs that can be formed from the current
     * list.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i=0; i<outerList.getAtomCount(); i++) {
            atoms.atom0 = outerList.getAtom(i);
            for (int j=0; j<innerList.getAtomCount(); j++) {
                atoms.atom1 = innerList.getAtom(j);
                action.actionPerformed(atoms);
            }
        }
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return outerList.getAtomCount() * innerList.getAtomCount();
    }

    /**
     * Returns 2, indicating that this is a pair iterator
     */
    public int nBody() {
        return 2;
    }

    /**
     * Sets the list that will be used to generate the pairs. Must call reset()
     * before beginning iteration.
     * 
     * @param atomList
     *            the new atom list for iteration
     */
    public void setOuterList(AtomSet newList) {
        this.outerList = newList;
        unset();
    }

    /**
     * Sets the list that will be used to generate the pairs. Must call reset()
     * before beginning iteration.
     * 
     * @param atomList
     *            the new atom list for iteration
     */
    public void setInnerList(AtomSet newList) {
        this.innerList = newList;
        unset();
    }

    /**
     * Returns the outer list used to generate the pairs.
     */
    public AtomSet getOuterList() {
        return outerList;
    }

    /**
     * Returns the inner list used to generate the pairs.
     */
    public AtomSet getInnerList() {
        return innerList;
    }

    private static final long serialVersionUID = 1L;
    private AtomSet outerList, innerList;
    private int outerIndex, innerIndex;
    private final AtomPair atoms = new AtomPair();
}
