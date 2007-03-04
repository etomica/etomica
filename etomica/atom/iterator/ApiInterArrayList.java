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
public class ApiInterArrayList implements AtomPairIterator, java.io.Serializable {

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
    public ApiInterArrayList(AtomArrayList outerList, AtomArrayList innerList) {
        if (outerList == innerList) {
            throw new IllegalArgumentException(
                    "ApiInterList will not work if both iterators are the same instance");
        }
        setOuterList(outerList);
        setInnerList(innerList);
    }

    /**
     * Returns true if the given pair of atoms are both in the current lists,
     * with the first atom considered from the outer list, and the second atom
     * from the inner list.
     */
    public boolean contains(AtomSet pair) {
        if (pair == null || pair.count() != 2 ||
                pair.getAtom(0) == pair.getAtom(1)) {
            return false;
        }
        return outerList.contains(pair.getAtom(0))
                && innerList.contains(pair.getAtom(1));
    }

    /**
     * Indicates whether iterator has another iterate.
     */
    public boolean hasNext() {
        return nextOuterIndex < outerList.size() && nextInnerIndex < innerList.size();
    }

    /**
     * Sets iterator in condition to begin iteration.
     * 
     * @throws IllegalStateException
     *             if outer and inner lists have been set to the same instance
     */
    public void reset() {
        if (outerList == innerList && outerList != emptyList) {
            throw new IllegalStateException(
                    "ApiInterList will not work correctly if inner and outer lists are the same instance");
        }
        nextOuterIndex = 0;
        nextInnerIndex = 0;
    }

    /**
     * Sets iterator such that hasNext is false.
     */
    public void unset() {
        nextOuterIndex = outerList.size();
    }

    /**
     * Same as nextPair.
     */
    public AtomSet next() {
        return nextPair();
    }

    /**
     * Returns the next iterate pair. Returns null if hasNext() is false.
     */
    public AtomPair nextPair() {
        if (!hasNext())
            return null;
        atoms.atom0 = outerList.get(nextOuterIndex);
        atoms.atom1 = innerList.get(nextInnerIndex);
        nextInnerIndex++;
        if (nextInnerIndex == innerList.size()) {
            nextOuterIndex++;
            nextInnerIndex = 0;
        }
        return atoms;
    }

    /**
     * Returns the next iterate pair without advancing the iterator.
     */
    public AtomSet peek() {
        atoms.atom0 = outerList.get(nextOuterIndex);
        atoms.atom1 = innerList.get(nextInnerIndex);
        return atoms;
    }

    /**
     * Performs given action on all pairs that can be formed from the current
     * list.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i=0; i<outerList.size(); i++) {
            atoms.atom0 = outerList.get(i);
            for (int j=0; j<innerList.size(); j++) {
                atoms.atom1 = innerList.get(j);
                action.actionPerformed(atoms);
            }
        }
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return outerList.size() * innerList.size();
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
    public void setOuterList(AtomArrayList newList) {
        this.outerList = (newList != null) ? newList : emptyList;
        unset();
    }

    /**
     * Sets the list that will be used to generate the pairs. Must call reset()
     * before beginning iteration.
     * 
     * @param atomList
     *            the new atom list for iteration
     */
    public void setInnerList(AtomArrayList newList) {
        this.innerList = (newList != null) ? newList : emptyList;
        unset();
    }

    /**
     * Returns the outer list used to generate the pairs.
     */
    public AtomArrayList getOuterList() {
        return outerList;
    }

    /**
     * Returns the inner list used to generate the pairs.
     */
    public AtomArrayList getInnerList() {
        return innerList;
    }

    private static final long serialVersionUID = 1L;
    private AtomArrayList outerList, innerList;
    private int nextOuterIndex, nextInnerIndex;
    private final AtomArrayList emptyList = new AtomArrayList();
    private final AtomPair atoms = new AtomPair();

}
