/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.action.MoleculeAction;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculeArrayList;
import etomica.atom.MoleculeSetSinglet;

/**
 * Class used to test iterators by collecting a list of strings
 * given by the toString method of the atoms (or atom sets) produced
 * by the iterator.
 * 
 * @author aawalker
 */
class ListerMolecule implements MoleculeAction {

	public final LinkedList<String> list;
	
	public ListerMolecule() {
		list = new LinkedList<String>();
	}
    
    public void clear() {
        list.clear();
    }

	public void actionPerformed(IMolecule atom) {
		list.add(atom.toString());
	}

	/**
     * Adds atomSet.toString() to list.
	 * kmb 4/27/05
	 */
	public void actionPerformed(IMoleculeList atomSet) {
        list.add(atomSet.toString());
	}
    
    /**
     * Performs action on each element of array.
     */
    public void addEachToList(IMoleculeList a) {
        for(int i=0; i<a.getMoleculeCount(); i++) {
            actionPerformed(new MoleculeSetSinglet(a.getMolecule(i)));
        }
    }

    /**
     * Performs action of each element of array.
     */
    public void addEachToList(MoleculeArrayList atomList) {
        MoleculeIteratorArrayListSimple iterator = new MoleculeIteratorArrayListSimple(atomList);
        iterator.reset();
        for (IMolecule atom = iterator.nextMolecule(); atom != null;
             atom = iterator.nextMolecule()) {
            actionPerformed(new MoleculeSetSinglet(atom));
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
    public static ListerMolecule[] listerArray(int n) {
        ListerMolecule[] lister = new ListerMolecule[n];
        for (int i=0;i<n;i++) {
            lister[i] = new ListerMolecule();
        }
        return lister;
    }


}
