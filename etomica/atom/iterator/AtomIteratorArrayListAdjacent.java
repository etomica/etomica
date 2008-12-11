package etomica.atom.iterator;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomToAtomLeafList;
import etomica.atom.AtomToIndex;
import etomica.atom.AtomToIndexChild;
import etomica.atom.AtomToParentChildList;
import etomica.atom.iterator.IteratorDirective.Direction;


/**
 * Returns one or both of the atoms adjacent to a specified atom in
 * its parent's child list. If adjacent linker has no atom no corresponding
 * atom is given as an iterate; thus iterator may give 0, 1, or 2 iterates
 * depending on presence of adjacent atoms and specification of iteration
 * direction.
 */
public class AtomIteratorArrayListAdjacent implements AtomIteratorAtomDependent, java.io.Serializable {

    /**
     * Constructor gives iterator not ready for iteration.  Must
     * set an atom and reset before using.  Default direction is
     * null, meaning that both adjacent atoms (if there are two)
     * will be given by iterator.
     */
    public AtomIteratorArrayListAdjacent(IteratorDirective.Direction direction) {
        this(direction,new AtomToIndexChild(), new AtomToParentChildList());
    }
    
    public AtomIteratorArrayListAdjacent(IteratorDirective.Direction direction,
            AtomToIndex atomToIndex, AtomToAtomLeafList atomToAtomSet) {
        super();
        this.direction = direction;
        this.atomToIndex = atomToIndex;
        this.atomToAtomSet = atomToAtomSet;
        unset();
        atomSetSinglet = new AtomSetSinglet();
    }

    /**
     * Returns the number of iterates that iterator would give
     * if reset and iterated in its current condition.  Does not
     * depend on or affect iteration state.
     */
    public int size() {
        int size = 0;
        if(direction != IteratorDirective.Direction.DOWN && (firstCursor < list.getAtomCount()-1)) {
            size++;
        }
        if(direction != IteratorDirective.Direction.UP && (firstCursor > 0)) {
            size++;
        }
        return size;
    }

    /**
     * Returns 1, indicating that this is an atom AtomSet iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Same as nextAtom.
     */
    public IAtomList next() {
        atomSetSinglet.atom = nextAtom();
        if (atomSetSinglet.atom == null) return null;
        return atomSetSinglet;
    }

    /**
     * Returns the next iterator, or null if hasNext is false.
     */
    public IAtomLeaf nextAtom() {
        if(upListNow) {
            upListNow = false;
            return list.getAtom(firstCursor+1);
        }
        else if (dnListNow) {
            dnListNow = false;
            return list.getAtom(firstCursor-1);
        }
        return null;
    }

    /**
     * Readies the iterator to begin iteration.
     */
    public void reset() {
        upListNow = ((direction != IteratorDirective.Direction.DOWN) && (firstCursor < list.getAtomCount()-1));
        dnListNow = ((direction != IteratorDirective.Direction.UP) && (firstCursor > 0));
    }

    /**
     * Puts iterator in a state where hasNext is false.
     */
    public void unset() {
        upListNow = false;
        dnListNow = false;
    }   
    
    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, depending on how iterator was configured at
     * construction.
     */
    public void setAtom(IAtomLeaf atom) {
        list = atomToAtomSet.getAtomList(atom);
        firstCursor = atomToIndex.getIndex(atom);
    }

    private static final long serialVersionUID = 1L;
    private int firstCursor;
    private final Direction direction;
    private boolean upListNow, dnListNow;
    private IAtomList list;
    private final AtomToIndex atomToIndex;
    private final AtomToAtomLeafList atomToAtomSet;
    protected final AtomSetSinglet atomSetSinglet;
}
