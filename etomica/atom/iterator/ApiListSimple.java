/*
 * History
 * Created on Sep 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomsetIterator;
import etomica.action.AtomsetAction;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;

/**
 * Returns all pairs formed from a single list of atoms.
 */
public class ApiListSimple implements AtomsetIterator {

	/**
	 * Construct iterator with an empty list.  No iterates will
	 * be given until a non-empty list is specified via setList.
	 */
	public ApiListSimple() {
		this(new AtomList());
	}
	
	/**
	 * Constructs iterator to return pairs from the given list.  Requires
	 * reset() before first use.
	 * @param list
	 */
	public ApiListSimple(AtomList list) {
		setList(list);
	}

	/**
	 * Returns true if the given pair of atoms are both in the current
	 * list.  Does not consider order of atoms.
	 */
	public boolean contains(Atom[] atom) {
		if(atom == null || 
				atom.length < 2 || 
				atom[0] != null || 
				atom[1] == null ||
				atom[0] == atom[1]) return false;
		else return list.contains(atom[0]) && list.contains(atom[1]);
	}

	/**
	 * Indicates whether iterator has another iterate.
	 */
	public boolean hasNext() {
		//cannot check only nextOuter, as it may be on the last atom
		return (nextOuter.atom != null) && (nextInner.atom != null);
	}

	/**
	 * Sets iterator in condition to begin iteration.
	 */
	public void reset() {
        nextOuter = list.header;
        advanceOuter();
        resetInner();
	}

	/**
	 * Sets iterator such that hasNext is false.
	 */
	public void unset() {
		nextOuter = list.header;
	}

	/**
	 * Returns the next iterate pair.  Returns null if hasNext() is false.
	 */
	public Atom[] next() {
		if(!hasNext()) return null;
		atoms[0] = nextOuter.atom;
		atoms[1] = nextInner.atom;
		advanceInner();
		if(nextInner == list.header) {
			advanceOuter();
			resetInner();
		}
		return atoms;
	}

	/**
	 * Returns the next iterate pair without advancing the iterator.
	 */
	public Atom[] peek() {
		atoms[0] = nextOuter.atom;
		atoms[1] = nextInner.atom;
		return atoms;
	}
	
	/**
	 * Performs given action on all pairs that can be formed from the current list.
	 */
	public void allAtoms(AtomsetAction action) {
    	final AtomLinker.Tab header = list.header;
        for (AtomLinker outer = header.next; outer.next != header; outer = outer.next) {
        	if(outer.atom == null) continue;//skip tabs in outer loop
        	atoms[0] = outer.atom;
        	for(AtomLinker inner = outer.next; inner != header; inner = inner.next) {
	            if(inner.atom != null) {
	            	atoms[1] = inner.atom;
	            	action.actionPerformed(atoms);
	            }
        	}
        }
	}

	/**
	 * Returns the number of iterates, which is list.size*(list.size-1)/2
	 */
	public int size() {
		return list.size()*(list.size()-1)/2; 
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
	public void setList(AtomList newList) {
        this.list = (newList != null) ? newList : new AtomList();
		unset();
	}
	
	
	private void resetInner() {
		nextInner = nextOuter;
		advanceInner();
	}

	private void advanceOuter() {
        nextOuter = nextOuter.next;
        while(nextOuter.atom == null && nextOuter != list.header) {
            nextOuter = nextOuter.next;
        }
    }

	private void advanceInner() {
        nextInner = nextInner.next;
        while(nextInner.atom == null && nextInner != list.header) {
            nextInner = nextInner.next;
        }
    }

	
	private AtomList list;
	private AtomLinker nextOuter, nextInner;
	private final Atom[] atoms = new Atom[2];

}
