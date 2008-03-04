package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomsetArrayList;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.species.Species;

/**
 * Iterator for all the molecules of a set of species in a box.  Each iterate
 * is all the molecules in a box, with each Atom as the first atom in the 
 * set. This class is used by PotentialMaster to iterate over molecules for 
 * N-body potentials.
 * 
 * This class is designed to work and conform to the API... not to be efficient 
 * or pleasant to look at!  Use neighbor lists. 
 */
public class AtomIteratorAllLeafType implements AtomsetIteratorPDT, java.io.Serializable {

    /**
     * @param species species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public AtomIteratorAllLeafType(IAtomType[] atomType) {
    	this.atomType = atomType;
        next = new AtomsetArrayList();
    }

    /**
     * Sets the box containing the molecules for iteration. A null
     * box conditions iterator to give no iterates.
     */
    public void setBox(IBox newBox) {
        box = newBox;
        if (box == null) {
            throw new NullPointerException("Null box");
        }
    }

    /**
     * Sets the target of iteration... has no actual effect since all iterates
     * contain all Atoms.
     */
    public void setTarget(IAtom newTargetAtom) {
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     * Besides, you didn't really want to iterate down, did you?
     */
    public void setDirection(Direction newDirection) {
    }

    public void reset() {
        // add all Atoms to ArrayList we will return
        AtomArrayList arrayList = next.getArrayList();
        arrayList.clear();
        IAtomSet leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
        	for (int j=0; j<atomType.length; j++) {
        		if(leafList.getAtom(i).getType()==atomType[j]){
        			arrayList.add(leafList.getAtom(i));
        		}
        	}
        }
        nextCursor = 0;
    }
    
    public void unset() {
        next.getArrayList().clear();
    }
    
    public IAtomSet next() {
        if (nextCursor + 1 > next.getAtomCount()) {
            return null;
        }
        if (nextCursor < 0) {
            // already poked
            nextCursor = -nextCursor;
            return next;
        }
        AtomArrayList arrayList = next.getArrayList();
        IAtom oldFirst = arrayList.getAtom(0);
        arrayList.set(0,arrayList.getAtom(nextCursor));
        arrayList.set(nextCursor,oldFirst);
        nextCursor++;
        return next;
    }
    
    public int nBody() {
        return Integer.MAX_VALUE;
    }
    
    public void allAtoms(AtomsetAction action) {
        reset();
        for (IAtomSet atoms = next(); atoms != null; atoms = next()) {
            action.actionPerformed(atoms);
        }
    }

    /**
     * Returns the number of iterates given by this iterator, if iterated after
     * a call to reset().
     */
    public int size() {
        return next.getAtomCount();
    }

    private static final long serialVersionUID = 1L;
    private final IAtomType[] atomType;
    private IBox box;
    private int nextCursor;
    private final AtomsetArrayList next;
}
