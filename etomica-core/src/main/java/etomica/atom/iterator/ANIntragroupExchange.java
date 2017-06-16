/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;

/**
 */
public class ANIntragroupExchange implements AtomsetIteratorBasisDependent {

    public ANIntragroupExchange(int nBody, int P) {
        atoms = new AtomArrayList(nBody);
        this.nBody = nBody;
        this.P = P;
        particleCounters = new int[nBody];
        done = new int[nBody];
        
        unset();
    }

    public void unset() {
        counter = -1;
    }

    public int size() {
    	if (nBody == 2) {
    		return P*particles/2*(particles-1); 
    	} else if (nBody == 3) {
    		return P*particles/3*(particles-2)/2*(particles-1); 
    	} else {
    		throw new RuntimeException("Not generalized for four-body nonadditive energies");
    	}
        
    }

    public final int nBody() {
        return nBody;
    }

    public IAtomList next() {
        if (done[0]==1) {
            return null;
        }
        atoms.clear();
        for (int i=0; i<nBody; i++) {
        	int particle = particleCounters[i]*P;
            atoms.add(atomList.getAtom(particle+counter));
        }
        //System.out.println(atoms);
        counter++;
        if (counter == P) {
        	counter=0;
        	for (int i=(nBody-1); i>=0; i--) {
        		if (particleCounters[i] > i+particles-nBody-1){
        			done[i]=1;
        			
        			
        		} else {
        			particleCounters[i]++;
        			for (int j=i+1; j<nBody; j++) {
        				particleCounters[j]=particleCounters[i]+1;
        			}
        			break;
        		}

        	}
        }
        return atoms;
    }

    public void setTarget(IAtom newTargetAtom) {
        // we are basis dependent because we have setBasis, but this
        // method should not be called.
    	if (newTargetAtom != null ){
    		throw new RuntimeException("nope");
    	}
        
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
    public void setBasis(IMoleculeList basisMolecule) {
        if (basisMolecule.getMoleculeCount() != 1) {
            throw new IllegalArgumentException();
        }
        atomList = basisMolecule.getMolecule(0).getChildList();
        particles = atomList.getAtomCount()/P;
        
    }

    public void reset() {
        counter = 0;
        for (int i=0;i<nBody;i++) {
        	particleCounters[i] = i;
        	done[i]=0;
        }
    }

    /**
     * Returns 2, indicating that the setBasis method expects an array of two
     * atoms.
     */
    public int basisSize() {
        return 1;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomArrayList atoms;
    protected int counter, particles;
    protected int[] particleCounters, done;
    protected IAtomList atomList;
    protected boolean needSetupIterators = true;
    protected final int nBody;
    protected final int P;
}
