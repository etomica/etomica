/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.atom.AtomPair;

/**
 * Iterator that expires after returning a single atom pair, which is specified
 * by a call to the setPair method. Subsequent calls to reset() and next() will
 * return the specified pair, until another is specified via setPair. No
 * iteration is performed if either or both atoms are null.
 * 
 * @author Tai Tan
 */
public class ApiIndexList implements AtomsetIteratorBasisDependent {

	/**
     * Constructs iterator without defining atoms in pair.
     */
    public ApiIndexList(int [][]index) {
        pair = new AtomPair();
        this.index = index;
    }

    public int basisSize(){
    	return 1;
    }
    
    public boolean haveTarget(IAtom a){
    	for(int i =0; i < index.length; i++){   //index.length = number of pairs
    		
    		if (a == parentGroup.getChildList().getAtom(index[i][0])){
    			return true;
    		}
    		if (a == parentGroup.getChildList().getAtom(index[i][1])){
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
        if (parentGroup != null) {
            cursor = 0;
        }
    }

    /**
     * Returns the iterator's pair and unsets iterator.
     */
    public IAtomList next() {
    	if (target != null){
        	for(; cursor < index.length; cursor++){   //index.length = number of pairs
        		
        		if (target == parentGroup.getChildList().getAtom(index[cursor][0])){
        			break;
        		}
        		if (target == parentGroup.getChildList().getAtom(index[cursor][1])){
        			break;
        		}
        	}
    	}
    	
        if (cursor >= index.length){
        	return null;
        }
		pair.atom0 = parentGroup.getChildList().getAtom(index[cursor][0]);
		pair.atom1 = parentGroup.getChildList().getAtom(index[cursor][1]);
		cursor++;
        return pair;
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public final int nBody() {
        return 2;
   }

    private int [][]index;
    private int cursor;
    private IMolecule parentGroup;
    private IAtom target;
    
    private final AtomPair pair;
    
	private static final long serialVersionUID = 1L;

}

