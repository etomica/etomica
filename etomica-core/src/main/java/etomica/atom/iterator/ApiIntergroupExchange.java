/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;

/**
 */
public class ApiIntergroupExchange extends ApiIntergroup {

    public ApiIntergroupExchange(int P) {
    	super(new AtomIteratorBasisAdjustableStart(), new AtomIteratorBasisAdjustableStart());
    	this.P = P;
    }


    /**
     * Returns the number of pairs given by this iterator. Not dependent on
     * state of hasNext.
     */
    public int size() {
    	int monomersOuter = aiOuter.size()/P;
    	int monomersInner = aiInner.size()/P;
        return P*monomersOuter*monomersInner;
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public IAtomList next() {
    	//Advance through first monomer of polymer1 and all of polymer 2; then second monomer of polymer1 and all of polymer 2, etc.
    
    	//outer = iterator over atoms in polymer 1 
    	//inner = iterator over atoms in polymer 2
    	
    	if (atom1 >= aiOuterSize) {
            return null;
        }
    	
        pair.atom0 = aiOuter.nextAtom();
        pair.atom1 = aiInner.nextAtom();
        //System.out.println(pair);
        atom1++;
        atom2++;
        if (atom2 == aiInnerSize) { 
        	atom2=0;
        	aiInner.reset();
        	monomer1++;
        	atom1=monomer1*P;
        	if (atom1 < aiOuterSize) {
        		((AtomIteratorAtomDependent)aiOuter).setAtom(((AtomIteratorArrayListSimple)aiOuter).getList().getAtom(atom1));
        	}
        	
        } else if (atom1 % P == 0) { 
        	atom1=monomer1*P;
        	((AtomIteratorAtomDependent)aiOuter).setAtom(((AtomIteratorArrayListSimple)aiOuter).getList().getAtom(atom1));
        	aiOuter.reset();
        }
        
        
        return pair;
    }

    public void reset() {
    	((AtomIteratorAtomDependent)aiOuter).setAtom(startAtom0);
    	((AtomIteratorAtomDependent)aiInner).setAtom(startAtom1);
    	aiOuterSize = aiOuter.size();
        aiInnerSize = aiInner.size();
        aiOuter.reset();
        aiInner.reset();
        monomer1 = 0;
        atom1 = 0;
        atom2 = 0;
    }
    
    public void setBasis(IMoleculeList basisAtoms) {
        super.setBasis(basisAtoms);
        if (basisAtoms != null) {
        	startAtom0 = basisAtoms.getMolecule(0).getChildList().getAtom(0);
        	startAtom1 = basisAtoms.getMolecule(1).getChildList().getAtom(0);
        }
    }
    
    protected final int P;
    protected int monomer1;
    protected int atom1;
    protected int atom2;
    protected IAtom startAtom0;
    protected IAtom startAtom1;
    protected int aiOuterSize;
    protected int aiInnerSize;
}
