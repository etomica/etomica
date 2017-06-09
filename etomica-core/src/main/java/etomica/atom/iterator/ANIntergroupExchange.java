/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMoleculeList;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculePair;
import etomica.atom.MoleculeSetSinglet;

/**
 */
public class ANIntergroupExchange implements AtomsetIteratorBasisDependent {

    public ANIntergroupExchange(int P) {
        atoms = new AtomArrayList(3);
        moleculePair = new MoleculePair();
        anIntragroupExchange = new ANIntragroupExchange(2,P);
        this.P = P;
        unset();
    }

    public void unset() {
        counter = -1;
    }

    public int size() {
    	int particles0 = moleculePair.atom0.getChildList().getAtomCount()/P;
    	int particles1 = moleculePair.atom1.getChildList().getAtomCount()/P;
    	
        return P*particles0/2*(particles0-1)*particles1+P*particles1/2*(particles1-1)*particles0;
    }

    public final int nBody() {
        return 3;
    }

    public IAtomList next() {
    	//find all corresponding triplets from two polymers: one from polymer 1 and two from polymer 2, and vice versa
    	// if possible, we begin by extracting pairs from polymer 0
    	atoms.clear();
    	
    	int size0 = moleculePair.atom0.getChildList().getAtomCount();
    	int size1 = moleculePair.atom1.getChildList().getAtomCount();
    	
    	
    	
    	if (stillOnPolymer0) {
    		if (counter == size1) {
    			return null;
    		}
    		atoms.addAll(pair);
    		atoms.add(moleculePair.atom1.getChildList().getAtom(counter));
    		counter++;
    		
    		if (counter%P==0) {
    			counter = counter - P; //counter tracking monomer on polymer 1
        	} 
    		//Have we exhausted the pairs on polymer 0?
    		pair = anIntragroupExchange.next();
    		if (pair == null) {
    			counter = counter + P; //counter tracking monomer on polymer 1
    			anIntragroupExchange.reset();
    			if (counter == size1) {
    				//switch to taking pairs from polymer1, if it is big enough
    				if (moleculePair.atom1.getChildList().getAtomCount()==P){
    					//all done with monomers from polymer 1
    					//System.out.println(atoms);
    					return atoms;
    				}
    				stillOnPolymer0 = false;
    				counter=0;
    			} 
    			
    		}
    		
    	} else {
    		if (counter == size0) {
    			return null;
    		}
    		atoms.add(moleculePair.atom0.getChildList().getAtom(counter));
    		atoms.addAll(pair);
    		counter++;
    		
    		if (counter%P==0) {
    			counter = counter - P; //counter tracking monomer on polymer 0
        	} 
    		//Have we exhausted the pairs on polymer 1?
    		pair = anIntragroupExchange.next();
    		if (pair == null) {
    			counter = counter + P; //counter tracking monomer on polymer 0
    			anIntragroupExchange.reset();
    			if (counter == size0) {
    				//all done
    				//System.out.println(atoms);
    				return atoms;
    			} 
    			
    		}
    	}
    	
    	
    	//System.out.println(atoms);
        
        return atoms;
    }

    public void setTarget(IAtom newTargetAtom) {
        // we are basis dependent because we have setBasis, but this
        // method should not be called.
        throw new RuntimeException("nope");
    }

    public boolean haveTarget(IAtom target) {
        // we are basis dependent because we have setBasis, but this
        // method should not be called.
        throw new RuntimeException("nope");
    }

    /**
     * Specifies the basis, which identifies the molecules subject to
     * iteration.
     */
    public void setBasis(IMoleculeList basisMolecules) {
        if (basisMolecules.getMoleculeCount() != 2) {
            throw new IllegalArgumentException();
        }
        
        moleculePair.atom0 = basisMolecules.getMolecule(0);
        moleculePair.atom1 = basisMolecules.getMolecule(1);
    }

    public void reset() {
        counter = 0;
        stillOnPolymer0 = moleculePair.atom0.getChildList().getAtomCount()>P;
        MoleculeSetSinglet one = new MoleculeSetSinglet(stillOnPolymer0? moleculePair.atom0 : moleculePair.atom1);
        anIntragroupExchange.setBasis(one);
        anIntragroupExchange.reset();
        pair = anIntragroupExchange.next();
    }

    /**
     * Returns 2, indicating that the setBasis method expects an array of two
     * atoms.
     */
    public int basisSize() {
        return 2;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomArrayList atoms;
    protected final MoleculePair moleculePair;
    protected int counter, nLeafAtoms;
    protected boolean needSetupIterators = true;
    protected final int P;
    protected ANIntragroupExchange anIntragroupExchange;
    protected boolean stillOnPolymer0;
    protected IAtomList pair;
}
