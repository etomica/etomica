package etomica.junit;

import etomica.*;
import java.util.LinkedList;

/**
 * @author aawalker
 *
 */
class Lister implements AtomsetActive, AtomActive {
	
	public final LinkedList list;
	
	public Lister() {
		list = new LinkedList();
	}

	/**
	 * @see etomica.AtomActive#actionPerformed(etomica.Atom)
	 */
	public void actionPerformed(Atom atom) {
		list.add(atom.toString());
	}

	/**
	 * @see etomica.AtomActive#actionPerformed(etomica.Atom)
	 */
	public void actionPerformed(Atom[] atom) {
		String newString = "[";
		for(int i=0; i<atom.length-1; i++) newString += atom[i].toString()+",";
		if(atom.length > 0) newString += atom[atom.length-1].toString();
		newString += "]";
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
