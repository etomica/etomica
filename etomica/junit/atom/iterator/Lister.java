package etomica.junit.atom.iterator;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.action.AtomAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorArrayListSimple;

/**
 * Class used to test iterators by collecting a list of strings
 * given by the toString method of the atoms (or atom sets) produced
 * by the iterator.
 * 
 * @author aawalker
 */
class Lister extends AtomsetActionAdapter implements AtomAction {

	public final LinkedList list;
	public AtomSet atoms;
	
	public Lister() {
		list = new LinkedList();
	}
    
    public void clear() {
        list.clear();
    }

	public void actionPerformed(IAtom atom) {
		list.add(atom.toString());
	}

	public void setAtom(IAtom atom) {
		setAtoms(atoms); //AtomSet={atoms});
	}

	public IAtom getAtom() {
		return (atoms == null || atoms.getAtomCount() < 1) ? null : atoms.getAtom(0);
	}

	/**
     * Adds atomSet.toString() to list.
	 * kmb 4/27/05
	 */
	public void actionPerformed(AtomSet atomSet) {
        list.add(atomSet.toString());
	}
    
    /**
     * Performs action on each element of array.
     */
    public void addEachToList(AtomSet a) {
        for(int i=0; i<a.getAtomCount(); i++) {
            actionPerformed(new AtomSetSinglet(a.getAtom(i)));
        }
    }

    /**
     * Performs action of each element of array.
     */
    public void addEachToList(AtomArrayList atomList) {
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(atomList);
        iterator.reset();
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            actionPerformed(new AtomSetSinglet(atom));
        }
    }
    
	public String toString() {
	    Iterator iter = list.iterator();
        String string = "{";
        while(iter.hasNext()) {
            string += iter.next().toString();
            if(iter.hasNext()) string += ",";
        }
        string += "}";
        return string;
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
