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
    
    public Class atomSequencerClass() {return IteratorFactorySimple.INSTANCE.atomSequencerClass();}
    
    public Class neighborSequencerClass() {return Sequencer.class;}
    
    public AtomSequencer.Factory atomSequencerFactory() {return IteratorFactorySimple.INSTANCE.atomSequencerFactory();}
    
    public AtomSequencer.Factory neighborSequencerFactory() {return Sequencer.FACTORY;}
    
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
    public boolean hasNext() {return listIterator.hasNext();}
    
    /**
     * True if the parent group of the given atom is the current basis for the iterator.
     * False otherwise, or if atom or basis is null.
     */
    public boolean contains(Atom atom) {
        return atom != null && basis != null && atom.node.parentNode() == basis;
    }
    
	public void setAsNeighbor(boolean b) {
	    throw new RuntimeException("method IteratorFactoryCell.IntraGroupIterator.setAsNeighbor not implemented");
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
    
    //we assume that the only Index links in the list are those demarking
    //the beginning of each cell's sequence; thus we reset the list iterator
    //using null as the terminator
    
    public Atom reset(Atom atom) {
        referenceAtom = atom;
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        nextAtom = null;
        if(atom == null) {
            throw new IllegalArgumentException("Cannot reset IteratorFactoryCell.IntragroupIterator without referencing an atom");
        //probably need isDescendedFrom instead of parentGroup here
        } 
        if(atom.node.parentNode() != basis) {
            throw new IllegalArgumentException("Cannot return IteratorFactoryCell.IntragroupIterator referencing an atom not in group of basis");
        }
        if(iterateCells) {
            referenceCell = (AtomCell)((Sequencer)atom.seq).site();
            listIterator.unset();
            if(upListNow) {
                cellIterator.reset(referenceCell, IteratorDirective.UP);//set cell iterator to return reference cell
                AtomCell cell = (AtomCell)cellIterator.next();//this returns referenceCell
                listIterator.reset(referenceAtom.seq, null, IteratorDirective.UP);
                listIterator.next();//advance so not to return reference atom
            }
            if(!listIterator.hasNext()) advanceCell();
        } else if(upListNow) {
            listIterator.reset(referenceAtom.seq, IteratorDirective.UP);
            listIterator.next();//advance so not to return reference atom
            if(!listIterator.hasNext() && doGoDown) {
                listIterator.reset(referenceAtom.seq, IteratorDirective.DOWN);
                listIterator.next();//advance so not to return reference atom
                upListNow = false;
                doGoDown = false;
            }
        } else if(doGoDown) {
            listIterator.reset(referenceAtom.seq, IteratorDirective.DOWN);
            listIterator.next();//advance so not to return reference atom
            doGoDown = false;
        }
        return listIterator.peek();
    }
                
    // Moves to next cell that has an iterate
    private void advanceCell() {
        do {
            if(cellIterator.hasNext()) {
                AtomCell cell = (AtomCell)cellIterator.next();
                if(upListNow) {
                    listIterator.reset(cell.first(speciesIndex), null, IteratorDirective.UP);
                } else {
                    listIterator.reset(cell.last(speciesIndex), null, IteratorDirective.DOWN);
                }
            } else if(doGoDown) {//no more cells that way; see if should now reset to look at down-cells
                cellIterator.reset(referenceCell, IteratorDirective.DOWN);//set cell iterator to return reference cell
                AtomCell cell = (AtomCell)cellIterator.next();//this will be the reference cell
                listIterator.reset(referenceAtom.seq, null, IteratorDirective.DOWN);
                listIterator.next();//advance so not to return reference atom
                upListNow = false;
                doGoDown = false;
            } else {//no more cells at all
                break;
            }
        } while(!listIterator.hasNext());
    }
            
    public Atom next() {
        Atom atom = listIterator.next();
        if(!listIterator.hasNext() && iterateCells) advanceCell();
        return atom;
    }
    
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
     * a leaf atom, a class-cast exception will be thrown.
     */
    public void setBasis(Atom atom) {
        setBasis((AtomTreeNodeGroup)atom.node);
    }
    
    public void setBasis(AtomTreeNodeGroup node) {
        basis = node;
        //can base decision whether to iterate over cells on type of sequencer
        //for given atom, because it is in the group of atoms being iterated
        iterateCells = basis.childSequencerClass().equals(Sequencer.class);
        if(iterateCells) speciesIndex = basis.parentSpecies().index;
        listIterator.setBasis(node.childList);
    }
    
    /**
     * Returns the current iteration basis.
     */
    public Atom getBasis() {return basis.atom();}
    
    /**
     * The number of atoms returned on a full iteration, using the current basis.
     */
    public int size() {return (basis != null) ? basis.childAtomCount() : 0;}   

    private AtomTreeNodeGroup basis;
    private Atom next;
    private Atom referenceAtom, nextAtom;
    private boolean upListNow, doGoDown;
    private IteratorDirective.Direction direction, currentDirection;
    private AtomCell referenceCell;
    private boolean iterateCells;
    private int speciesIndex;
    private final SiteIteratorNeighbor cellIterator = new SiteIteratorNeighbor();
    private final AtomIteratorList listIterator = new AtomIteratorList();

}//end of IntragroupIterator
   
/////////////////////////////////////////////////////////////////////////////////////////////

public static final class Sequencer extends AtomSequencer implements AbstractLattice.Occupant {
    
    public AtomCell cell;             //cell currently occupied by this coordinate
    public Lattice lattice;           //cell lattice in the phase occupied by this coordinate
    
    public Sequencer(Atom a) {
        super(a);
    }

    public Site site() {return cell;}   //Lattice.Occupant interface method

    public void setLattice(Lattice newLattice) {
        lattice = newLattice;
        if(lattice != null) assignCell();
    }

    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
     //this methods needs to be fixed
    public boolean preceeds(Atom a) {
        //want to return false if atoms are the same atoms
        if(a == null) return true;
        if(atom.node.parentGroup() == a.node.parentGroup()) {
            if(((Sequencer)atom.seq).site().equals(cell)) {
                //this isn't correct
                return atom.node.index() < a.node.index();//works also if both parentGroups are null
            }
            else return ((Sequencer)atom.seq).site().preceeds(cell);
        }
        int thisDepth = atom.node.depth();
        int atomDepth = a.node.depth();
        if(thisDepth == atomDepth) return atom.node.parentGroup().seq.preceeds(a.node.parentGroup());
        else if(thisDepth < atomDepth) return this.preceeds(a.node.parentGroup());
        else /*if(this.depth > atom.depth)*/ return atom.node.parentGroup().seq.preceeds(a);
    }

//Determines appropriate cell and assigns it
    public void assignCell() {
        AtomCell newCell = (AtomCell)lattice.nearestSite(atom.coord.position());
        if(newCell != cell) {assignCell(newCell);}
    }
//Assigns atom to given cell; if removed from another cell, repairs tear in list
    public void assignCell(AtomCell newCell) {
        cell = newCell;
        if(cell == null) return;
        AtomTreeNodeGroup parentNode = atom.node.parentNode();
        parentNode.childList.moveBefore(this, newCell.last(parentNode.parentSpecies().index));
    }//end of assignCell
    
    public static final AtomSequencer.Factory FACTORY = new AtomSequencer.Factory() {
        public AtomSequencer makeSequencer(Atom atom) {return new Sequencer(atom);}
    };
}//end of Sequencer

/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * A factory that makes Sites of type AtomCell
 */
private static final class AtomCellFactory implements SiteFactory {
    public Site makeSite(AbstractLattice parent, AbstractLattice.Coordinate coord) {
        if(!(coord instanceof BravaisLattice.Coordinate)) {
            throw new IllegalArgumentException("IteratorFactoryCell.AtomCellFactory: coordinate must be of type BravaisLattice.Coordinate");
        }
        AbstractCell cell = ((Lattice)parent).makeWignerSeitzCell((AbstractLattice.PositionCoordinate)coord);
        return new AtomCell((Lattice)parent, cell);
    }
}//end of AtomCellFactory
    
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * A lattice cell that holds a reference to a sequence of atoms.
 */
private static final class AtomCell extends AbstractCell {
    public Space.Vector position;
    AbstractCell cell;
    private AtomLinker.Index[] firstIndex, lastIndex;
    public AtomCell(Lattice parent, AbstractCell cell) {
        super(parent, (AbstractLattice.PositionCoordinate)cell.coordinate());
        this.cell = cell;
//        color = Constants.RandomColor();
//            position = (Space2D.Vector)coord.position();
    }
    public boolean inCell(Space.Vector r) {return cell.inCell(r);}
    public int D() {return cell.D();}
    public Space.Vector[] vertex() {return cell.vertex();}
    public double volume() {return cell.volume();}
    public AtomLinker.Index first(int speciesIndex) {return firstIndex[speciesIndex];}
    public AtomLinker.Index last(int speciesIndex) {return lastIndex[speciesIndex];}
}//end of AtomCell

   
}//end of IteratorFactoryCell