package etomica.atom.iterator;

import etomica.api.IAtomList;

/**
 * Iterator that returns coupled iterates; the first pair contains the first
 * atom from each basis.  The second pair has the second atom from each basis.
 *
 * @author Andrew Schultz
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
    	
    	if (atom1 >= aiOuter.size()) {
            return null;
        }
    	
        pair.atom0 = aiOuter.nextAtom();
        pair.atom1 = aiInner.nextAtom();
        atom1++;
        atom2++;
        if (atom2 == aiInner.size()) { 
        	atom2=0;
        	aiInner.reset();
        	monomer1++;
        	atom1=monomer1*P;
        	if (atom1 < aiOuter.size()) {
        		((AtomIteratorAtomDependent)aiOuter).setAtom(((AtomIteratorArrayListSimple)aiOuter).getList().getAtom(atom1));
        	}
        	
        } else if (atom1 % P == 0) { 
        	atom1=monomer1*P;
        	((AtomIteratorAtomDependent)aiOuter).setAtom(((AtomIteratorArrayListSimple)aiOuter).getList().getAtom(atom1));
        }
        
        
        return pair;
    }

    public void reset() {
        aiOuter.reset();
        aiInner.reset();
        monomer1 = 0;
        atom1 = 0;
        atom2 = 0;
    }
    
    protected final int P;
    protected int monomer1;
    protected int atom1;
    protected int atom2;
}
