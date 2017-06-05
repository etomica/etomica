/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * General class for a collection of linearly linked atoms.
 * 
 * @author nancycribbin
 */

public abstract class ConformationChain implements IConformation, java.io.Serializable {

    public ConformationChain(Space space){
		this.space = space;
		//orientationVector = space.makeVector();
		//wrongNumberOfVectors = "Wrong number of vectors in the argument to ConformationChain subclass.";
	}
	
	/**
	 * Resets the vector instructions.
	 */
	protected abstract void reset();
	
	/**
	 * 
	 * @return the instructions to get to the location of the next molecule from
	 * the current one.
	 */
	protected abstract Vector nextVector();
	
	/**
	 * Places a set of atoms in a linearly connected fashion.
	 * @param atomlist a list of atoms in the order in which they are linked
	 */
	public void initializePositions(IAtomList atomlist){
		
		//First, check that we actually have some atoms
		int size = atomlist.getAtomCount();
    	if(size == 0) return;
    
    	reset();
	
    	//space.makeVector() zeroes the made Vector automatically
    	Vector currentPosition = space.makeVector();
    
    	//Zero the first atom.
        atomlist.getAtom(0).getPosition().E(0);
    	
        for (int iLeaf=1; iLeaf<size; iLeaf++) {
            IAtom a = atomlist.getAtom(iLeaf);
    		//TODO someday, we might want a to be a chunk-of-atoms
    		currentPosition.PE(nextVector());
    		a.getPosition().E(currentPosition);
    	}
	}
	
    private static final long serialVersionUID = 1L;
    protected final Space space;
	/**
	 * The vector drawn from the head of the molecule to the tail of the molecule.
	 */
/*	private Vector calculateOrientationVector(AtomList atomlist){
		//atomIterator.reset();
		Atom a = atomlist.getFirst();
		Atom z = atomlist.getLast();
		
		Vector ov = space.makeVector();
		
		ov.E(z.coord.position());
		ov.ME(a.coord.position());
		return ov;
	}
*/	
	
	//TODO Should we have a method here that moves the atom somehow- rotate, translate, etc.?

	//private Vector orientationVector;	
}
