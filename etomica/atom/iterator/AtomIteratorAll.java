package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetDetect;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomsetArrayList;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Iterator for all the molecules of a set of species in a phase.  Each iterate
 * is all the molecules in a phase, with each Atom as the first atom in the 
 * set. This class is used by PotentialMaster to iterate over molecules for 
 * N-body potentials.
 * 
 * This class is designed to work and conform to the API... not to be efficient 
 * or pleasant to look at!  Use neighbor lists. 
 */
public class AtomIteratorAll implements AtomsetIteratorPDT, java.io.Serializable {

    /**
     * @param species species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public AtomIteratorAll(Species[] species) {
        this.species = species;
        next = new AtomsetArrayList();
    }

    /**
     * Sets the phase containing the molecules for iteration. A null
     * phase conditions iterator to give no iterates.
     */
    public void setPhase(Phase newPhase) {
        phase = newPhase;
    }

    /**
     * Sets the target of iteration... has no actual effect since all iterates
     * contain all Atoms.
     */
    public void setTarget(Atom newTargetAtom) {
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
        for (int i=0; i<species.length; i++) {
            AtomTreeNodeGroup speciesAgentNode = (AtomTreeNodeGroup)phase.getAgent(species[i]).node;
            arrayList.addAll(speciesAgentNode.childList);
        }
    }
    
    public boolean hasNext() {
        return nextCursor < next.count();
    }
    
    public void unset() {
        next.getArrayList().clear();
    }
    
    public AtomSet next() {
        if (!hasNext()) {
            return null;
        }
        if (nextCursor < 0) {
            // already poked
            nextCursor = -nextCursor;
            return next;
        }
        AtomArrayList arrayList = next.getArrayList();
        Atom oldFirst = arrayList.get(0);
        arrayList.set(0,arrayList.get(nextCursor));
        arrayList.set(nextCursor,oldFirst);
        nextCursor++;
        return next;
    }
    
    public int nBody() {
        return Integer.MAX_VALUE;
    }
    
    public void allAtoms(AtomsetAction action) {
        reset();
        while (hasNext()) {
            action.setAtoms(next());
            action.actionPerformed();
        }
    }

    public AtomSet peek() {
        if (!hasNext()) {
            return null;
        }
        if (nextCursor < 0) {
            // already poked
            return next;
        }
        AtomArrayList arrayList = next.getArrayList();
        Atom oldFirst = arrayList.get(0);
        arrayList.set(0,arrayList.get(nextCursor));
        arrayList.set(nextCursor,oldFirst);
        // leave a secret Masonic signal to tell nextPair it doesn't need to 
        // advance
        nextCursor = -(nextCursor + 1);
        return next;
    }

    /**
     * Returns the number of iterates given by this iterator, if iterated after
     * a call to reset().
     */
    public int size() {
        return next.count();
    }

    /**
     * Returns true if given atom set is among the iterates, if iterated after a
     * call to reset().
     */
    public boolean contains(AtomSet atom) {
        if (atom == null)
            return false;
        detector.setAtoms(atom);
        detector.reset();
        allAtoms(detector);
        return detector.detectedAtom();
    }

    private final AtomsetDetect detector = new AtomsetDetect(null);
    private final Species[] species;
    private Phase phase;
    private int nextCursor;
    private final AtomsetArrayList next;
}
