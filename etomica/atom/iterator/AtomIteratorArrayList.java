package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
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
 	protected Atom atoms[] = new Atom[1];
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
 	
 	public Atom[] next() {
 		atoms[0] = nextAtom();
 		return atoms;
 	}
 
 	public Atom[] peek() {
 		atoms[0] = list.get(cursor);
 		return atoms;
 	}
 	public int size() {
 		return list.size();
 	}
 	
 	public void allAtoms(AtomsetAction act) {
 		int arraySize = size();
 		for (int i=0; i<arraySize; i++) {
 			atoms[0] = list.get(i);
 			act.actionPerformed(atoms);
 		}
 	}
 	
 	public void reset() {
 		cursor = 0;
 	}
 	
 	public boolean contains(Atom[] atom) {
 		return contains(atom);
 	}
 
 }