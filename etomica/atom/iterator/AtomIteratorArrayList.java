package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.atom.AtomArrayList;

 /**
  * An atom iterator of the elements from an AtomArrayList (in proper
  * sequence).
  */

public class AtomIteratorArrayList implements AtomIterator {
	/**
 	 * Index of element to be returned by subsequent call to next.
 	 */
 	protected int cursor = 0;
 	protected AtomArrayList list;

 	public AtomIteratorArrayList() {
 		list = new AtomArrayList();
 	}
 	public AtomIteratorArrayList(AtomArrayList atomList) {
 		list = atomList;
 	}
 	public void setList(AtomArrayList atomList) {
 		list = atomList;
 	}
 	
 	public int nBody() {return 1;}
 	public void unset() {cursor = size();}
 
 	public boolean hasNext() {
 	    return cursor != size();
 	}
 
 	public Atom nextAtom() {
 		return list.get(cursor++);
 	}
 	
 	public AtomSet next() {
 		return nextAtom();
 	}
 
 	public AtomSet peek() {
 		return list.get(cursor);
 	}
 	public int size() {
 		return list.size();
 	}
 	
 	public void allAtoms(AtomsetAction act) {
 		int arraySize = size();
 		for (int i=0; i<arraySize; i++) {
 			act.actionPerformed(list.get(i));
 		}
 	}
 	
 	public void reset() {
 		cursor = 0;
 	}
 	
 	public boolean contains(AtomSet atom) {
 		return contains(atom);
 	}
 
 }