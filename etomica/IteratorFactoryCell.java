package etomica;

import etomica.lattice.*;

public class IteratorFactoryCell implements IteratorFactory {
    
    public static final IteratorFactoryCell INSTANCE = new IteratorFactoryCell();
    
    public AtomIterator makeAtomIterator() {return new AtomIteratorChildren();}
        
    public AtomIterator makeIntragroupIterator() {return new IntragroupIterator();}
    public AtomIterator makeIntergroupIterator() {return new AtomIteratorChildren();}
    
    public AtomSequencer makeAtomSequencer(Atom atom) {
        return IteratorFactorySimple.INSTANCE.makeAtomSequencer(atom);
    }
    public AtomSequencer makeNeighborSequencer(Atom atom) {return new Sequencer(atom);}
    //maybe need an "AboveNbrLayerSequencer" and "BelowNbrLayerSequencer"
    
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Iterates among the children of a given basis, those atoms
 * that are cell-list neighbors of a specified atom that is
 * a child of the same basis.
 */
//would like to modify so that central atom can be any descendant of the basis.
public static final class IntragroupIterator implements AtomIterator {
    
    /**
     * Indicates if another iterate is forthcoming.
     */
    public boolean hasNext() {return nextAtom != null;}
    
    /**
     * True if the parent group of the given atom is the current basis for the iterator.
     * False otherwise, or if atom or basis is null.
     */
    public boolean contains(Atom atom) {
        return atom != null && basis != null && atom.node.parentGroup() == basis;
    }
    
    /**
     * Does reset if atom in iterator directive is child of the current basis.  
     * Sets hasNext false if given atom does is not child of basis.  Throws
     * an IllegalArgumentException if directive does not specify an atom.
     */
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        return reset(id.atom1());
    }
    
    public Atom reset(Atom atom) {
        referenceAtom = atom;
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        nextAtom = null;
        if(atom == null) {
            throw new IllegalArgumentException("Cannot reset IteratorFactoryCell.IntragroupIterator without referencing an atom");
        //probably need isDescendedFrom instead of parentGroup here
        } 
        if(atom.node.parentGroup() != basis) {
            throw new IllegalArgumentException("Cannot return IteratorFactoryCell.IntragroupIterator referencing an atom not in grouup of basis");
        }
        //can base decision whether to iterate over cells on type of sequencer
        //for given atom, because it is in the group of atoms being iterated
        iterateCells = (atom.seq instanceof Sequencer);
        if(iterateCells) {
            cell = ((Sequencer)atom.seq).site();
            if(upListNow) {
                cellIterator.reset(cell, IteratorDirective.UP);//set cell iterator to return next cell up (shouldn't begin with this cell)
                nextAtom = atom.seq.nextAtom();
            }
            if(nextAtom == null) advanceCell();
        } else if(upListNow) {
            nextAtom = atom.seq.nextAtom();
        } else if(doGoDown) {
            nextAtom = atom.seq.previousAtom();
        }
        return nextAtom;
    }
                
    // Finds first atom of next occupied cell
    private void advanceCell() {
        do {
            if(cellIterator.hasNext() && iterateCells) {
                nextAtom = upListNow ? ((AtomCell)cellIterator.next()).first()
                                     : ((AtomCell)cellIterator.next()).last();
            } else if(doGoDown) {//no more cells that way; see if should now reset to look at down-cells
                cellIterator.reset(cell, IteratorDirective.DOWN);//set cell iterator to return next cell down
                nextAtom = referenceAtom.seq.previousAtom();
                upListNow = false;
                doGoDown = false;
            } else {//no more cells at all
                break;
            }
        } while(nextAtom == null);
    }
            
    public Atom next() {
        Atom atom = nextAtom;
        nextAtom = upListNow ? atom.seq.nextAtom() : atom.seq.previousAtom();
        if(nextAtom == null) advanceCell();
        return atom;
    }
    /**
     * Ignored.
     */
    public void setAsNeighbor(boolean b) {}
    
    /**
     * Throws RuntimeException because this is a neighbor iterator, and must
     * be reset with reference to an atom.
     */
    public Atom reset() {
        throw new RuntimeException("Cannot reset IteratorFactoryCell.IntragroupIterator without referencing an atom");
    }
    
    
    /**
     * Performs given action for each child atom of basis.
     */
    public void allAtoms(AtomAction act) {
        throw new RuntimeException("AtomIteratorNbrCellIntra.allAtoms not implemented");
/*        if(basis == null) return;
        last = basis.node.lastChildAtom();
        for(Atom atom = basis.node.firstChildAtom(); atom != null; atom=atom.nextAtom()) {
            act.actionPerformed(atom);
            if(atom == last) break;
        }*/
    }
        
    /**
     * Sets the given atom as the basis, so that child atoms of the
     * given atom will be returned upon iteration.  If given atom is
     * a leaf atom, no iterates are given.
     */
    public void setBasis(Atom atom) {
        basis = /*(atom == null || atom.node.isLeaf()) ? null :*/ atom;
    }
    
    /**
     * Returns the current iteration basis.
     */
    public Atom getBasis() {return basis;}
    
    /**
     * The number of atoms returned on a full iteration, using the current basis.
     */
    public int size() {return (basis != null) ? basis.node.childAtomCount() : 0;}   

    private Atom basis;
    private Atom next;
    private Atom referenceAtom;
    private boolean upListNow, doGoDown;
    private IteratorDirective.Direction direction;
    private AtomCell cell;
    private boolean iterateCells;

}//end of IntragroupIterator
   
/////////////////////////////////////////////////////////////////////////////////////////////

public static final class Sequencer implements AtomSequencer, AbstractLattice.Occupant {
    
    Atom next, previous;  //next and previous coordinates in neighbor list
    public AtomCell cell;             //cell currently occupied by this coordinate
    public Lattice lattice;               //cell lattice in the phase occupied by this coordinate
    private final Atom atom;
    
    public Sequencer(Atom a) {
        atom = a;
    }
    public void setNextAtom(Atom a) {
        next = a;
        if(a != null) a.seq.setPreviousAtom(this.atom);
    }
    public void setPreviousAtom(Atom a) {
        previous = a;
    }
    public Site site() {return cell;}   //Lattice.Occupant interface method
    public void clearPreviousAtom() {previous = null;}
    public Atom nextAtom() {return next;}
    public Atom previousAtom() {return previous;}

    public void setLattice(Lattice newLattice) {
        lattice = newLattice;
        if(lattice != null) assignCell();
    }

//Determines appropriate cell and assigns it
    public void assignCell() {
        AtomCell newCell = (AtomCell)lattice.nearestSite(atom.position(), dimensions);
        if(newCell != cell) {assignCell(newCell);}
    }
//Assigns atom to given cell; if removed from another cell, repairs tear in list
    public void assignCell(AtomCell newCell) {
        if(previous != null) {previous.setNextAtom(next);}
        else {//removing first atom in cell
            if(cell != null) cell.setFirst(next); 
            if(next != null) next.clearPreviousAtom();
        }   
        cell = newCell;
        if(cell == null) {setNextAtom(null);}
        else {
            setNextAtom(cell.first);
            cell.first = this;
        }
        clearPreviousAtom();
    }//end of assignCell
}//end of Sequencer

/**
 * A factory that makes Sites of type AtomCell
 */
private static class AtomCellFactory implements SiteFactory {
    public Site makeSite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
        if(!(coord instanceof BravaisLattice.Coordinate)) {
            throw new IllegalArgumentException("IteratorFactoryCell.AtomCellFactory: coordinate must be of type BravaisLattice.Coordinate");
        }
        return (new AtomCell(parent, (BravaisLattice.Coordinate)coord));
    }
}//end of AtomCellFactory
    
/**
 * A lattice cell that holds a reference to an atom coordinate, which marks the beginning of list of atoms in that cell.
 */
private static class AtomCell extends SquareLattice.Cell {
    public Coordinate first;
    public Space2D.Vector position;
    public Color color;
    public AtomCell N, E, S, W;   //immediately adjacent cells (North, East, South, West)
    public double yN, xE, yS, xW; //x, y coordinates of boundaries of cell
    public AtomCell(Lattice parent, BravaisLattice.Coordinate coord) {
        super(parent, coord);
        color = Constants.RandomColor();
//            position = (Space2D.Vector)coord.position();
    }
    public Coordinate first() {return first;}
    public void setFirst(Coordinate f) {first = f;}
}//end of AtomCell

   
}//end of IteratorFactoryCell