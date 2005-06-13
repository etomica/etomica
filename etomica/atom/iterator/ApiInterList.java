/*
 * History
 * Created on Sep 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;

/**
 * Returns all pairs formed from two different lists of atoms.
 * 
 */
public class ApiInterList implements AtomPairIterator {

	/**
	 * Construct iterator with an empty lists.  No iterates will
	 * be given until non-empty lists are specified via setList.
	 */
	public ApiInterList() {
		this(new AtomList(), new AtomList());
	}
	
	/**
	 * Constructs iterator to return pairs from the given list.  Requires
	 * reset() before first use.
	 * @param list
	 */
	public ApiInterList(AtomList outerList, AtomList innerList) {
		setOuterList(outerList);
        setInnerList(innerList);
	}

	/**
	 * Returns true if the given pair of atoms are both in the current
	 * list.  Does not consider order of atoms.
	 */
	public boolean contains(AtomSet pair) {
		if(!(pair instanceof AtomPair) || 
				((AtomPair)pair).atom0 == ((AtomPair)pair).atom1) return false;
		return outerList.contains(((AtomPair)pair).atom0) 
                                && innerList.contains(((AtomPair)pair).atom1);
	}

	/**
	 * Indicates whether iterator has another iterate.
	 */
	public boolean hasNext() {
		return (nextOuter.atom != null) && (nextInner.atom != null);
	}

	/**
	 * Sets iterator in condition to begin iteration.
	 */
	public void reset() {
        nextOuter = outerList.header.next;
        nextInner = innerList.header.next;
	}

	/**
	 * Sets iterator such that hasNext is false.
	 */
	public void unset() {
		nextOuter = outerList.header;
	}

    public AtomSet next() {
        return nextPair();
    }
    
	/**
	 * Returns the next iterate pair.  Returns null if hasNext() is false.
	 */
	public AtomPair nextPair() {
		if(!hasNext()) return null;
		atoms.atom0 = nextOuter.atom;
		atoms.atom1 = nextInner.atom;
		nextInner = nextInner.next;
		if(nextInner.atom == null) {
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
	 * Performs given action on all pairs that can be formed from the current list.
	 */
	public void allAtoms(AtomsetAction action) {
        for (AtomLinker outer = outerList.header.next; outer.atom != null; outer = outer.next) {
        	atoms.atom0 = outer.atom;
        	for(AtomLinker inner = innerList.header.next; inner.atom != null; inner = inner.next) {
            	atoms.atom1 = inner.atom;
            	action.actionPerformed(atoms);
        	}
        }
	}

	/**
	 * Returns the number of iterates, which is list.size*(list.size-1)/2
	 */
	public int size() {
		return outerList.size()*innerList.size();
	}

	/**
	 * Returns 2, indicating that this is a pair iterator
	 */
	public int nBody() {
		return 2;
	}
	
	/**
	 * Sets the list that will be used to generate the pairs.
	 * Must call reset() before beginning iteration.
	 * @param atomList the new atom list for iteration
	 */
	public void setOuterList(AtomList newList) {
        this.outerList = (newList != null) ? newList : emptyList;
		unset();
	}
    
    /**
     * Sets the list that will be used to generate the pairs.
     * Must call reset() before beginning iteration.
     * @param atomList the new atom list for iteration
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
