package etomica.junit;

import java.util.LinkedList;

import etomica.Atom;
import etomica.action.AtomAction;
import etomica.action.AtomsetAction;

/**
 * @author aawalker
 *
 */
class Lister implements AtomsetAction, AtomAction {
	
	public final LinkedList list;
	
	public Lister() {
		list = new LinkedList();
	}

	/**
	 * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
	 */
	public void actionPerformed(Atom atom) {
		list.add(atom.toString());
	}

	/**
	 * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
	 */
	public void actionPerformed(Atom[] atom) {
//		for(int i=0; i<atom.length; i++) actionPerformed(atom[i]);
//		String newString = "[";
		String newString = "";
		for(int i=0; i<atom.length-1; i++) newString += atom[i].toString()+",";
		if(atom.length > 0) newString += atom[atom.length-1].toString();
//		newString += "]";
		list.add(newString);
	}

	/**
	 * Gives an array of listers each with their own list.
	 * 
	 * @param n the number of listers in the array
	 * @return AtomLister[]
	 */	
	public static Lister[] listerArray(int n) {
		Lister[] lister = new Lister[n];
		for (int i=0;i<n;i++) {
			lister[i] = new Lister();
		}
		return lister;
	}

}
