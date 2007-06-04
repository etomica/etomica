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
 * Returns all pairs formed from a single list of atoms.
 */
public class ApiIntraArrayList implements AtomsetIterator, java.io.Serializable {

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
     * Sets iterator in condition to begin iteration.
     */
    public void reset() {
        if (list.getAtomCount() < 2) {
            outerIndex = 2;
            innerIndex = 2;
            return;
        }
        outerIndex = 0;
        innerIndex = 0;
        atoms.atom0 = list.getAtom(0);
    }

    /**
     * Sets iterator such that next is null.
     */
    public void unset() {
        outerIndex = list.getAtomCount() - 2;
        innerIndex = list.getAtomCount() - 1;
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return list.getAtomCount() * (list.getAtomCount() - 1) / 2;
    }

    /**
     * Returns the next iterate pair. Returns null if hasNext() is false.
     */
    public AtomSet next() {
        if (innerIndex > list.getAtomCount() - 2) {
            if (outerIndex > list.getAtomCount() - 3) {
                return null;
            }
            outerIndex++;
            atoms.atom0 = list.getAtom(outerIndex);
            innerIndex = outerIndex;
        }
        innerIndex++;
        atoms.atom1 = list.getAtom(innerIndex);
        return atoms;
    }

    /**
     * Performs given action on all pairs that can be formed from the current
     * list.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i=0; i<list.getAtomCount()-1; i++) {
            atoms.atom0 = list.getAtom(i);
            for (int j=i+1; j<list.getAtomCount(); j++) {
                atoms.atom1 = list.getAtom(j);
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
        list = newList;
        unset();
    }

    /**
     * Returns the list used to generate the pairs.
     */
    public AtomArrayList getList() {
        return list;
    }

    private static final long serialVersionUID = 1L;
    private AtomArrayList list;
    private int outerIndex, innerIndex;
    private final AtomPair atoms = new AtomPair();
}
