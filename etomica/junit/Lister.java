package etomica.junit;

import etomica.*;
import etomica.action.*;
import java.util.LinkedList;

/**
 * @author aawalker
 *
 */
class Lister extends AtomsetActionAdapter implements AtomAction {
// Below is original line, from kmb2 project as of 12/22/04 work
//class Lister implements AtomsetActive, AtomActive {	

	public final LinkedList list;
	public AtomSet atoms;
	
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
/*	public void actionPerformed(Atom[] atom) {
//		for(int i=0; i<atom.length; i++) actionPerformed(atom[i]);
//		String newString = "[";
		String newString = "";
		System.out.println("atom.length is: "+atom.length); // added 1/7/05 to debug Atom
		for(int i=0; i<atom.length-1; i++) newString += atom[i].toString()+",";
		if(atom.length > 0) newString += atom[atom.length-1].toString();
//		newString += "]";
		list.add(newString);
	}
*/
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

	public void setAtom(Atom atom) {
//		setAtoms(new Atom[] {atom}); // could use this structure to wrap the atom during actionperformed; use with arrays only, not single atom in IteratorTest?
		setAtoms(atoms); //AtomSet={atoms});
		}

	public Atom getAtom() {
		return (atoms == null || atoms.count() < 1) ? null : atoms.getAtom(0);
	}

	/* New actionPerformed method based on AtomSet rather than Atom[] as argument.
	 * kmb 4/27/05
	 */
	public void actionPerformed(AtomSet a) {
		String newString = "";
		System.out.println("atom.length is: " + a.count());
		for(int i=0; i<a.count()-1; i++) newString += a.getAtom(i).toString()+",";
		if(a.count() > 0) newString += a.getAtom(a.count()-1).toString();

		list.add(newString);
	}



}
