package etomica.atom.iterator;


import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.atom.AtomsetArray;

/**
 * Atomset Iterator that iterates over set-of-four atoms
 * 
 * @author Tai Tan
 */


public class Atomset4IteratorIndexList implements AtomsetIteratorBasisDependent, AtomsetIterator {

	
	
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
                                

	public void setBasis(IAtomSet parent) {
	    if (parent == null) {
	        parentGroup = null;
	    }
	    else {
	        parentGroup = (IMolecule)parent.getAtom(0);
	    }
		unset();
	}

    public int size() {
        return index.length;
    }

    /**
     * Performs the given action on an AtomPair containing the three atoms last
     * identified via setSet. Does nothing if either such atom is null.
     * Unaffected by and has no effect on the reset/unset state.
     */
    public void allAtoms(AtomsetAction action) {
        if (parentGroup == null) {
            return;
        }
        
    	for(int i =0; i < index.length; i++){   //index.length = number of sets
    		
    		atoms[0] = parentGroup.getChildList().getAtom(index[i][0]);
    		atoms[1] = parentGroup.getChildList().getAtom(index[i][1]);
    		atoms[2] = parentGroup.getChildList().getAtom(index[i][2]);
    		atoms[3] = parentGroup.getChildList().getAtom(index[i][3]);
    		action.actionPerformed(atomset);
    	}
    }

    /**
     * Returns true if three non-null atoms have set and a call to reset() has
     * been performed, without any subsequent calls to next() or nextPair().
     */
    public boolean hasNext() {

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

    /**
     * Returns the iterator's pair and unsets iterator.
     */
    public AtomsetArray nextSet() {
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
     * Same as nextSet().
     */
    public IAtomSet next() {
        return nextSet();
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

