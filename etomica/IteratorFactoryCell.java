package etomica;
import etomica.utility.HashMap;

import etomica.lattice.*;

public class IteratorFactoryCell implements IteratorFactory {
    
    private Primitive primitive;
    private Simulation simulation;
    private int[] dimensions;
    private BravaisLattice[] deployedLattices = new BravaisLattice[0];
    
    public IteratorFactoryCell(Simulation sim) {
        this(sim, new PrimitiveCubic(sim.space), 2);
    }
    
    public IteratorFactoryCell(Simulation sim, Primitive primitive, int nCells) {
        this.simulation = sim;
        this.primitive = primitive;
        dimensions = new int[sim.space.D()];
        for(int i=0; i<sim.space.D(); i++) dimensions[i] = nCells;
    }
    
    public BravaisLattice makeCellLattice(final Phase phase) {
        AtomFactory cellFactory = primitive.unitCellFactory();
        ((PrimitiveCubic)primitive).setSize(phase.boundary().dimensions().component(0)/(double)dimensions[0]);//this needs work
        AtomFactory latticeFactory = new BravaisLattice.Factory(simulation.space, cellFactory, dimensions, primitive);
        BravaisLattice lattice = (BravaisLattice)latticeFactory.build();
        NeighborManager.Criterion neighborCriterion = new NeighborManager.Criterion() {
            public boolean areNeighbors(Site s1, Site s2) {
                return ((AbstractCell)s1).r2NearestVertex((AbstractCell)s2, phase.boundary()) < 0.2;
            }
        };
        lattice.setupNeighbors(neighborCriterion);
        lattice.agents = new Object[1];
        lattice.agents[0] = new HashMap();
        
        //resize and update lattice array as needed
        int latticeCountOld = deployedLattices.length;
        if(phase.index >= latticeCountOld) {
            BravaisLattice[] newArray = new BravaisLattice[phase.index+1];
            for(int i=0; i<latticeCountOld; i++) newArray[i] = deployedLattices[i];
            deployedLattices = newArray;
        }
        deployedLattices[phase.index] = lattice;
        
        return lattice;
    }
    
    public BravaisLattice getLattice(Phase phase) {
        if(phase == null) return null;
        else return deployedLattices[phase.index];
    }
    
    public AtomIterator makeGroupIteratorSimple() {return new FixedSequenceIterator();}

    public AtomIterator makeGroupIteratorSequential() {return new AtomIteratorChildren();}
        
    public AtomIterator makeIntragroupIterator() {return new IntragroupIterator(this);}
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
    
    public IntragroupIterator(IteratorFactoryCell factory) {
        iteratorFactory = factory;
    }
    
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
    
    //we assume that the only Tab links in the list are those demarking
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
            referenceCell = (AbstractCell)((Sequencer)atom.seq).site();
            cellIterator.setBasis(referenceCell);
            listIterator.unset();
            if(upListNow) {
                cellIterator.reset(IteratorDirective.UP);//set cell iterator to return first up-neighbor of reference cell
                listIterator.reset(referenceAtom.seq, null, IteratorDirective.UP);
                listIterator.next();//advance so not to return reference atom
            }
            if(!listIterator.hasNext()) advanceCell();
        } else if(upListNow) {//no cell iteration
            listIterator.reset(referenceAtom.seq, IteratorDirective.UP);
            listIterator.next();//advance so not to return reference atom
            if(!listIterator.hasNext() && doGoDown) {
                listIterator.reset(referenceAtom.seq, IteratorDirective.DOWN);
                listIterator.next();//advance so not to return reference atom
                upListNow = false;
                doGoDown = false;
            }
        } else if(doGoDown) {//no cell iteration
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
                Atom cell = cellIterator.next();
                AtomLinker.Tab[] tabs = (AtomLinker.Tab[])cell.agents[0];
                if(upListNow) {
                    listIterator.reset(tabs[speciesIndex], null, IteratorDirective.UP);
                } else {
                    listIterator.reset(tabs[speciesIndex].nextTab, null, IteratorDirective.DOWN);
                }
            } else if(doGoDown) {//no more cells that way; see if should now reset to look at down-cells
                cellIterator.reset(IteratorDirective.DOWN);//set cell iterator to return first down neighbor of reference cell
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
    
    //may be room for efficiency here
    public void setBasis(AtomTreeNodeGroup node) {
        basis = node;
        //can base decision whether to iterate over cells on type of sequencer
        //for given atom, because it is in the group of atoms being iterated
        iterateCells = basis.childSequencerClass().equals(Sequencer.class);
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
    private AbstractCell referenceCell;
    private boolean iterateCells;
    private int speciesIndex;
    private final SiteIteratorNeighbor cellIterator = new SiteIteratorNeighbor();
    private final AtomIteratorList listIterator = new AtomIteratorList();
    private final IteratorFactoryCell iteratorFactory;

}//end of IntragroupIterator
   
/////////////////////////////////////////////////////////////////////////////////////////////

public final class FixedSequenceIterator implements AtomIterator {
    
    private AtomList list;
	private Sequencer next;
	private AtomLinker first;
        
	public FixedSequenceIterator() {
	    this(AtomList.NULL);
	}
	public FixedSequenceIterator(AtomList list) {
	    setBasis(list);
	    reset();
	}
	
    //when iterator expires next will be set to null
	public boolean hasNext() {return next != null;}
	
	/**
	 * Sets the childList of the given atom as the basis for iteration.
	 * If atom is a leaf, hasNext is false.
	 */
    public void setBasis(Atom atom){
        if(atom.node.isLeaf()) unset();
        else setBasis(((AtomTreeNodeGroup)atom.node).childList);
    }
    public void setBasis(AtomList list) {
        if(list == null) list = AtomList.NULL;
        this.list = list;
        next = null;
    }
    public Atom getBasis(){return null;}//no implementation

    /**
     * Performs action on all atoms.
     */
    public void allAtoms(AtomAction action){
        reset();
        while(next != null) {
            action.actionPerformed(next.atom);
            next = next.nextFixed;
            if(next == first) next = null;
        }
    }
    
    /**
     * Sets iterator so that it is ready to go up its entire list of iterates.
     */
    public Atom reset() {
        first = list.header.next;
        while(first.atom == null && first != list.header) first = first.next;
        if(first == list.header) {next = null; return null;}
        else {next = (Sequencer)first; return next.atom;}
    }
    
    /**
     * Resets, ignoring given directive (same as reset()).
     */
    public Atom reset(IteratorDirective id) {return reset();}
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {next = null;}
    
    /**
     * Returns true if the given atom is in the list of iterates, false otherwise.
     */
	public boolean contains(Atom atom){
        return list.contains(atom);
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size() { return list.size();}

	    
    public Atom next() {
        return nextLinker().atom;
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.
     */
    public Atom peek() {
        return (next != null) ? next.atom : null;
    }
    
    private AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = next.nextFixed;
        if(next == first) next = null;
        return nextLinker;
    }//end of nextLinker

}//end of FixedSequenceIterator


/////////////////////////////////////////////////////////////////////////////////////////////

public static final class Sequencer extends AtomSequencer implements AbstractLattice.Occupant {
    
    public AbstractCell cell;                 //cell currently occupied by this coordinate
    public BravaisLattice lattice;    //cell lattice in the phase occupied by this coordinate
    private int listIndex;
    public Sequencer nextFixed, previousFixed;
    
    public Sequencer(Atom a) {
        super(a);
        nextFixed = previousFixed = this;
    }

    public Site site() {return cell;}   //Lattice.Occupant interface method

    public int listIndex() {return listIndex;}

    public void remove() {
        super.remove();
	    previousFixed.nextFixed = nextFixed;
	    nextFixed.previousFixed = previousFixed;
    }
        
    public void addBefore(AtomLinker newNext) {
        super.addBefore(newNext);
        while(newNext.atom == null && newNext != this) newNext = newNext.next;
        if(newNext == this) {//this is the first and only non-tab entry in the list
            nextFixed = previousFixed = this;
            return;
        }
        nextFixed = (Sequencer)newNext;
        previousFixed = nextFixed.previousFixed;
        previousFixed.nextFixed = this;
        nextFixed.previousFixed = this;
	}
	/**
	 * Reshuffles position of "standard" links without altering the fixed links.
	 */
	public void moveBefore(AtomLinker newNext) {
	    super.remove();
	    super.addBefore(newNext);
	}


    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
     //this method needs to be fixed
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
    
    /**
     * Method called when a translate method of coordinate is invoked.
     */
    public void moveNotify() {
        if(!cell.inCell(atom.coord.position())) assignCell();
    }
    
    /**
     * Method called when the parent of the atom is changed.
     * By the time this method is called, the atom has been placed
     * in the childList of the given parent (if it is not null).
     */
    public void setParentNotify(AtomTreeNodeGroup newParent) {
        if(newParent == null) {
            cell = null;
            lattice = null;
            return;
        }
        //get cell lattice for the phase containing the parent
        lattice = ((IteratorFactoryCell)newParent.parentSimulation().iteratorFactory).getLattice(newParent.parentPhase());
        //determine the index used by the cells for their tabs in the parent's childList
        HashMap hash = (HashMap)lattice.agents[0];
        Integer index = (Integer)hash.get(newParent);
        if(index == null) {//parent's childList isn't yet represented by cells
            index = new Integer(hash.size());
            hash.put(newParent, index);
            setupListTabs(lattice, newParent.childList);
        }
        listIndex = index.intValue();
        assignCell();
    }

//Determines appropriate cell and assigns it
    public void assignCell() {
        AbstractCell newCell = (AbstractCell)lattice.site(lattice.getPrimitive().latticeIndex(atom.coord.position()));
        if(newCell != cell) {assignCell(newCell);}
    }
//Assigns atom to given cell
    public void assignCell(AbstractCell newCell) {
        cell = newCell;
        if(cell == null) return;
        AtomTreeNodeGroup parentNode = atom.node.parentNode();
        parentNode.childList.moveBefore(this, ((AtomLinker.Tab[])newCell.agents[0])[listIndex].nextTab);
    }//end of assignCell
    
    public static final AtomSequencer.Factory FACTORY = new AtomSequencer.Factory() {
        public AtomSequencer makeSequencer(Atom atom) {return new Sequencer(atom);}
    };
}//end of Sequencer

/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * A factory that makes Sites of type AtomCell
 */
/*private static final class AtomCellFactory extends AtomFactory {
    public AtomCellFactory(Space space) {
        super(space);
    }
    public Atom build() {
        return new AtomCell();
    }
}//end of AtomCellFactory
    
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Wraps a lattice cell and holds a reference to a sequence of atoms.
 */
/*private static final class AtomCell extends AbstractCell {
    public final AbstractCell cell;
    private AtomLinker.Tab[] firstTab, lastTab;
    public AtomCell() {
//        this.cell = cell;
//        color = Constants.RandomColor();
//            position = (Space2D.Vector)coord.position();
    }
    public AtomLinker.Tab first(int speciesIndex) {return firstTab[speciesIndex];}
    public AtomLinker.Tab last(int speciesIndex) {return lastTab[speciesIndex];}
}//end of AtomCell
*/

/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Iterates through the cells of the given lattice, and adds a tab to the 
 * given list for each cell, and extends the tablist in each cell to reference
 * its new tab.
 */
private static void setupListTabs(BravaisLattice lattice, AtomList list) {
    AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
    iterator.reset();
    while(iterator.hasNext()) {
        Site site = (Site)iterator.next();
        if(site.agents == null || site.agents[0] == null) {
            site.agents = new Object[1];
            site.agents[0] = new AtomLinker.Tab[0];
        }
        AtomLinker.Tab[] tabList = (AtomLinker.Tab[])site.agents[0];
        AtomLinker.Tab[] newTabList = new AtomLinker.Tab[tabList.length+1];
        for(int i=0; i<tabList.length; i++) newTabList[i] = tabList[i];
        AtomLinker.Tab newTab = new AtomLinker.Tab();
        newTabList[tabList.length] = newTab;
        list.add(newTab);
        site.agents[0] = newTabList;
    }
}

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default.ATOM_SIZE = 1.0;
        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space2D());
        Simulation.instance = sim;
        IteratorFactoryCell iteratorFactory = new IteratorFactoryCell(sim);
        sim.setIteratorFactory(iteratorFactory);
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesSpheresMono speciesSpheres = new SpeciesSpheresMono();
	    speciesSpheres.setNMolecules(200);
	    Phase phase = new Phase();
	    iteratorFactory.makeCellLattice(phase);
	    Potential2 potential = new P2HardSphere();
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
        integratorHard.setTimeStep(0.01);
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
        
     //   controller.start();
    }//end of main
    
    
   
}//end of IteratorFactoryCell