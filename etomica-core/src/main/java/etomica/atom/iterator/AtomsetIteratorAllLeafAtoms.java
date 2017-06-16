/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomListWrapper;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;

/**
 * Iterator for all the molecules of a set of species in a box.  Each iterate
 * is all the molecules in a box, with each Atom as the first atom in the 
 * set. This class is used by PotentialMaster to iterate over molecules for 
 * N-body potentials.
 * 
 * This class is designed to work and conform to the API... not to be efficient 
 * or pleasant to look at!  Use neighbor lists. 
 */
public class AtomsetIteratorAllLeafAtoms implements AtomsetIteratorBasisDependent, java.io.Serializable {

    public AtomsetIteratorAllLeafAtoms() {
        next = new AtomListWrapper();
    }

    /**
     * Sets the target of iteration... has no actual effect since all iterates
     * contain all Atoms.
     */
    public void setTarget(IAtom newTargetAtom) {
    }

    public void reset() {
        AtomArrayList atomList = next.getArrayList();
    	atomList.clear();
    }
    
    public void unset() {
        next.getArrayList().clear();
    }
    
    public IAtomList next() {
    	if(next.getArrayList().getAtomCount()>0){
    		return null;
    	}
        for (int i=0; i<basis.getMoleculeCount(); i++){
        	next.getArrayList().add(basis.getMolecule(i).getChildList().getAtom(0));
        }
        return next;
    }
    
    public int nBody() {
        return Integer.MAX_VALUE;
    }
    
    /**
     * Returns the number of iterates given by this iterator, if iterated after
     * a call to reset().
     */
    public int size() {
        return next.getAtomCount();
    }

    private static final long serialVersionUID = 1L;
    private final AtomListWrapper next;
	private IMoleculeList basis;
	public int basisSize() {
		return Integer.MAX_VALUE;
	}

	public boolean haveTarget(IAtom target) {
		for (int i=0; i<basis.getMoleculeCount(); i++){
			if(basis.getMolecule(i).getChildList().getAtom(0) == target){
				return true;
			}
        }
		return false;
	}

	public void setBasis(IMoleculeList atoms) {
		basis = atoms;		
	}
}
