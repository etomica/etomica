package etomica.paracetamol;


import etomica.action.AtomsetAction;
import etomica.atom.AtomGroup;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;

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
                                

	public void setBasis(AtomSet parent) {
	    if (parent == null) {
	        parentGroup = null;
	    }
	    else {
	        parentGroup = (AtomGroup)parent.getAtom(0);
	    }
		unset();
	}

    public int size() {
        return index.length;
    }

    /**
     * Performs the given action on an AtomPair containing the two atoms last
     * identified via setPair. Does nothing if either such atom is null.
     * Unaffected by and has no effect on the reset/unset state.
     */
    public void allAtoms(AtomsetAction action) {
    	
    	for(int i =0; i < index.length; i++){   //index.length = number of pairs
    		pair.atom0 = parentGroup.getChildList().getAtom(index[i][0]);
    		pair.atom1 = parentGroup.getChildList().getAtom(index[i][1]);
    		
    		action.actionPerformed(pair);
    	}
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
    public AtomSet next() {
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
    private IAtomGroup parentGroup;
    private IAtom target;
    
    private final AtomPair pair;
    
	private static final long serialVersionUID = 1L;

}

