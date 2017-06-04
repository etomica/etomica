/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.atom.AtomsetArray;

/**
 * Atomset Iterator that iterates over set-of-four atoms
 * 
 * @author Tai Tan
 */


public class Atomset4IteratorIndexList implements AtomsetIteratorBasisDependent {
	
	/**
     * Constructs iterator without defining set of atoms.
     */
    public Atomset4IteratorIndexList(int [][]index) {
    	atomset = new AtomsetArray(4);
    	atoms = atomset.getArray();
        this.index = index;
    }

    public int basisSize(){
    	return 1;
    }
    
    public boolean haveTarget(IAtom a) {
        if (parentGroup == null) {
            return false;
        }
        
    	for(int i = 0; i < index.length; i++){   //index.length = number of sets
    		
    		if (a == parentGroup.getChildList().getAtom(index[i][0])){
    			return true;
    		}
    		if (a == parentGroup.getChildList().getAtom(index[i][1])){
    			return true;
    		}
    		if (a == parentGroup.getChildList().getAtom(index[i][2])){
    			return true;
    		}
    		if (a == parentGroup.getChildList().getAtom(index[i][3])){
    			return true;
    		}
    	}
    	return false;
    }
    
    public void setTarget(IAtom a){
    	target = a;
    	unset();
    }
                                

	public void setBasis(IMoleculeList parent) {
	    if (parent == null) {
	        parentGroup = null;
	    }
	    else {
	        parentGroup = parent.getMolecule(0);
	    }
		unset();
	}

    public int size() {
        return index.length;
    }

    /**
     * Returns true if three non-null atoms have set and a call to reset() has
     * been performed, without any subsequent calls to next() or nextPair().
     */
    protected boolean hasNext() {

    	if (target != null){
        	for(; cursor < index.length; cursor++){   //index.length = number of pairs
        		
        		if (target == parentGroup.getChildList().getAtom(index[cursor][0])){
        			break;
        		}
        		if (target == parentGroup.getChildList().getAtom(index[cursor][1])){
        			break;
        		}
        		if (target == parentGroup.getChildList().getAtom(index[cursor][2])){
        			break;
        		}
        		if (target == parentGroup.getChildList().getAtom(index[cursor][3])){
        			break;
        		}
        	}
    	}
        return (cursor < index.length);
    	
    }

    /**
     * Sets iterator to a state where hasNext() returns false.
     */
    public void unset() {
        cursor = index.length;
    }

    /**
     * Resets iterator to a state where hasNext is true, if atoms in pair are
     * not null.
     */
    public void reset() {
        if (parentGroup == null) {
            return;
        }
        
    	cursor = 0;
    }

    public IAtomList next() {
        if (!hasNext())
            return null;
        atoms[0] = parentGroup.getChildList().getAtom(index[cursor][0]);
        atoms[1] = parentGroup.getChildList().getAtom(index[cursor][1]);
        atoms[2] = parentGroup.getChildList().getAtom(index[cursor][2]);
        atoms[3] = parentGroup.getChildList().getAtom(index[cursor][3]);
        
        cursor++;
        return atomset;
    }

    /**
     * Returns 4, indicating that this is a set of four iterator.
     */
    public final int nBody() {
        return 4;
   }

    private int [][]index;
    private int cursor;
    private IMolecule parentGroup;
    private IAtom target;
    
    private AtomsetArray atomset;
    private IAtom[] atoms;
    
	private static final long serialVersionUID = 1L;

}

