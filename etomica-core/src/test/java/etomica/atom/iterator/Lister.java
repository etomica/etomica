/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.action.AtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSetSinglet;

/**
 * Class used to test iterators by collecting a list of strings
 * given by the toString method of the atoms (or atom sets) produced
 * by the iterator.
 * 
 * @author aawalker
 */
class Lister implements AtomAction {

	public final LinkedList<String> list;
	
	public Lister() {
		list = new LinkedList<String>();
	}
    
    public void clear() {
        list.clear();
    }

	public void actionPerformed(IAtom atom) {
		list.add(atom.toString());
	}

	/**
     * Adds atomSet.toString() to list.
	 * kmb 4/27/05
	 */
	public void actionPerformed(IAtomList atomSet) {
        list.add(atomSet.toString());
	}
    
    /**
     * Performs action on each element of array.
     */
    public void addEachToList(IAtomList a) {
        for(int i = 0; i<a.size(); i++) {
            actionPerformed(new AtomSetSinglet(a.get(i)));
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
	    Iterator<String> iter = list.iterator();
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
