/*
 * History
 * Created on Sep 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Returns all pairs formed from a single list of atoms. Can work with tabbed
 * list.
 */
public class ApiIntraArrayList implements AtomPairIterator, java.io.Serializable {

    /**
     * Construct iterator with an empty list. No iterates will be given until a
     * non-empty list is specified via setList.
     */
    public ApiIntraArrayList() {
        this(new AtomArrayList());
    }

    /**
     * Constructs iterator to return pairs from the given list. Requires reset()
     * before first use.
     * 
     * @param list
     */
    public ApiIntraArrayList(AtomArrayList list) {
        setList(list);
    }

    /**
     * Returns true if the given pair of atoms are both in the current list.
     * Does not consider order of atoms.
     */
    public boolean contains(AtomSet pair) {
        if (pair == null || pair.count() != 2 
                || pair.getAtom(0) == pair.getAtom(1))
            return false;
        return list.contains(pair.getAtom(0))
                && list.contains(pair.getAtom(1));
    }

    /**
     * Indicates whether iterator has another iterate.
     */
    public boolean hasNext() {
        //cannot check only nextOuter, as it may be on the last atom
        return (nextOuterIndex != -1) && (nextInnerIndex != -1);
    }

    /**
     * Sets iterator in condition to begin iteration.
     */
    public void reset() {
        nextOuterIndex = -1;
        advanceOuter();
        resetInner();
    }

    /**
     * Sets iterator such that hasNext is false.
     */
    public void unset() {
        nextOuterIndex = -1;
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return list.size() * (list.size() - 1) / 2;
    }

    public AtomSet next() {
        return nextPair();
    }

    /**
     * Returns the next iterate pair. Returns null if hasNext() is false.
     */
    public AtomPair nextPair() {
        if (!hasNext())
            return null;
        atoms.atom0 = list.get(nextOuterIndex);
        atoms.atom1 = list.get(nextInnerIndex);
        advanceInner();
        if (nextInnerIndex == -1) {
            advanceOuter();
            resetInner();
        }
        return atoms;
    }

    /**
     * Returns the next iterate pair without advancing the iterator.
     */
    public AtomSet peek() {
        atoms.atom0 = list.get(nextOuterIndex);
        atoms.atom1 = list.get(nextInnerIndex);
        return atoms;
    }

    /**
     * Performs given action on all pairs that can be formed from the current
     * list.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i=0; i<list.size()-1; i++) {
            atoms.atom0 = list.get(i);
            for (int j=i+1; j<list.size(); j++) {
                atoms.atom1 = list.get(j);
                action.actionPerformed(atoms);
            }
        }
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
    public void setList(AtomArrayList newList) {
        this.list = (newList != null) ? newList : new AtomArrayList();
        unset();
    }

    /**
     * Returns the list used to generate the pairs.
     */
    public AtomArrayList getList() {
        return list;
    }

    private void resetInner() {
        nextInnerIndex = nextOuterIndex;
        advanceInner();
    }

    private void advanceOuter() {
        nextOuterIndex++;
        if (nextOuterIndex == list.size()) {
            nextOuterIndex = -1;
        }
    }

    private void advanceInner() {
        nextInnerIndex++;
        if (nextInnerIndex == list.size()) {
            nextInnerIndex = -1;
        }
    }

    private AtomArrayList list;
    private int nextOuterIndex, nextInnerIndex;
    private final AtomPair atoms = new AtomPair();

}
