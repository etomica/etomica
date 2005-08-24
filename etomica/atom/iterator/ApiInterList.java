/*
 * History
 * Created on Sep 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Returns all pairs formed from two different untabbed lists of atoms.
 * Incorrect behavior will result if both lists refer to the same instance.
 *  
 */
public class ApiInterList implements AtomPairIterator, java.io.Serializable {

    /**
     * Construct iterator with an empty lists. No iterates will be given until
     * non-empty lists are specified via setList.
     */
    public ApiInterList() {
        this(new AtomList(), new AtomList());
    }

    /**
     * Constructs iterator to return pairs from the given lists. Requires
     * reset() before first use.
     * 
     * @throws IllegalArgumentException
     *             if both lists refer to the same instance
     */
    public ApiInterList(AtomList outerList, AtomList innerList) {
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
        return (nextOuter.atom != null) && (nextInner.atom != null);
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
        nextOuter = outerList.header.next;
        nextInner = innerList.header.next;
    }

    /**
     * Sets iterator such that hasNext is false.
     */
    public void unset() {
        nextOuter = outerList.header;
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
        atoms.atom0 = nextOuter.atom;
        atoms.atom1 = nextInner.atom;
        nextInner = nextInner.next;
        if (nextInner.atom == null) {
            nextOuter = nextOuter.next;
            nextInner = innerList.header.next;
        }
        return atoms;
    }

    /**
     * Returns the next iterate pair without advancing the iterator.
     */
    public AtomSet peek() {
        atoms.atom0 = nextOuter.atom;
        atoms.atom1 = nextInner.atom;
        return atoms;
    }

    /**
     * Performs given action on all pairs that can be formed from the current
     * list.
     */
    public void allAtoms(AtomsetAction action) {
        for (AtomLinker outer = outerList.header.next; outer.atom != null; outer = outer.next) {
            atoms.atom0 = outer.atom;
            for (AtomLinker inner = innerList.header.next; inner.atom != null; inner = inner.next) {
                atoms.atom1 = inner.atom;
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
    public void setOuterList(AtomList newList) {
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
    public void setInnerList(AtomList newList) {
        this.innerList = (newList != null) ? newList : emptyList;
        unset();
    }

    /**
     * Returns the outer list used to generate the pairs.
     */
    public AtomList getOuterList() {
        return outerList;
    }

    /**
     * Returns the inner list used to generate the pairs.
     */
    public AtomList getInnerList() {
        return innerList;
    }

    private AtomList outerList, innerList;
    private AtomLinker nextOuter, nextInner;
    private final AtomList emptyList = new AtomList();
    private final AtomPair atoms = new AtomPair();

}
