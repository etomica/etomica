package etomica;

import etomica.lattice.*;

/**
 * Variant of IteratorFactoryCell class, modified to work with leaf atoms
 * in particular, treating all equally in neighborlist.  Ignores species, basis,
 * and other features of atom hierarchy.
 * @author David Kofke
 * 03.08.10
 * 
 */

/* History
 * 08/10/03 (DAK) new, built from IteratorFactoryCell
 * 08/12/03 (DAK) modified as indicated in comments in code
 * 
 */
public class IteratorFactoryCellLeaf extends IteratorFactoryCell {
        
	/**
	 * Constructs a new iterator factory for the given simulation, using
	 * cubic cells for the neighbor listing.  Does not automatically
	 * register the factory with the simulation; this must be done
	 * separately using the simulation's setIteratorFactory method.
	 * For Simulation sim, at typical call is
	 * sim.setIteratorFactory(new IteratorFactoryCell(sim));
	 */
	public IteratorFactoryCellLeaf(Simulation sim) {
		super(sim);
	}
    
	/**
	 *  
	 * @param sim  The simulation in which this factory is being used
	 * @param nCells the number of cell in each dimension; total number of cells
	 * is then nCells^D.
	 */
	public IteratorFactoryCellLeaf(Simulation sim, int nCells) {
		super(sim, nCells);
	}
    
	/**
	 * Constructs a new iterator factory for the given simulation, using
	 * cells based on the given primitive.  Does not automatically
	 * register the factory with the simulation; this must be done
	 * separately using the simulation's setIteratorFactory method.
	 *
	 * @param sim          The simulation in which this factory is being used
	 * @param primitive    The primitive class that defines the type of unit
	 * cell used to construct the neighbor-cell lattice
	 * @param nCells the number of cell in each dimension; total number of cells is then nCells^D.
	 */
	public IteratorFactoryCellLeaf(Simulation sim, Primitive primitive, int nCells) {
		super(sim, primitive, nCells);
	}
    
	// use makeCellLattice method inherited from IteratorFactoryCell, even though lattice.agents is not used
	public BravaisLattice makeCellLattice(final Phase phase) {
		BravaisLattice lattice = super.makeCellLattice(phase);
		setupListTabs(lattice);
		return lattice;
	}

	public AtomIterator makeGroupIteratorSequential() {
		return new SequentialIterator(this);
	}
        
	public AtomIterator makeIntragroupNbrIterator() {return new IntragroupNbrIterator(this);}
	public AtomIterator makeIntergroupNbrIterator() {
		return null;//throw new RuntimeException("Inappropriate attempt to construct IntergroupNbrIterator; neighbor iterators loop over all leaf atoms, so there are no separate groups for iteration");
	}
    
	public AtomSequencer makeNeighborSequencer(Atom atom) {return new NeighborSequencer(atom);}
    
	public Class simpleSequencerClass() {return SimpleSequencer.class;}
    
	public Class neighborSequencerClass() {return NeighborSequencer.class;}
    
	public AtomSequencer.Factory simpleSequencerFactory() {return SimpleSequencer.FACTORY;}
    
	public AtomSequencer.Factory neighborSequencerFactory() {return NeighborSequencer.FACTORY;}
    
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Iterates all atoms, returning them in an order consistent with neighborlist
 * sequence but without limiting them to those neighboring a reference atom.
 */
public static final class SequentialIterator extends AtomIterator {
    
	private AtomIteratorList listIterator = new AtomIteratorList();
	private AtomList neighborSequenceList = new AtomList();
	private IteratorFactoryCell factory;
	private AtomTreeNodeGroup basis;
	private Phase lastPhase = null;
    
	public SequentialIterator(IteratorFactoryCell f) {
		factory = f;
	}
        
	public void all(Atom atom, IteratorDirective id, final AtomActive action) {
		if(atom == null || atom.node.isLeaf() || action == null) return;
		basis = (atom != null) ? (AtomTreeNodeGroup)atom.node : null;
		Phase phase = basis.parentPhase();
		BravaisLattice lattice = factory.getLattice(phase);
		Atom cell = lattice.siteList().getFirst();
		AtomLinker.Tab cellHeader = (AtomLinker.Tab)cell.agents[0];
		neighborSequenceList.setAsHeader(cellHeader, phase.speciesMaster.atomList.size());
		AtomIteratorList.all(cellHeader, id, action);		
	}

   /**
	 * Defines the atoms that are subject to iteration as the children of the
	 * given atom.
	 */
	public void setBasis(Atom a) {
		basis = (a != null) ? (AtomTreeNodeGroup)a.node : null;
		if(basis == null /*|| basis.childAtomCount() == 0*/) {
			listIterator.setBasis(AtomList.NULL);
			return;
		}
		boolean iterateCells = true;//basis.childSequencerClass().equals(NeighborSequencer.class);
		if(iterateCells) {
			Phase phase = basis.parentPhase();
			if(phase != lastPhase) {
				lastPhase = phase;
				BravaisLattice lattice = factory.getLattice(phase);
				Atom cell = lattice.siteList().getFirst();
				AtomLinker.Tab cellHeader = (AtomLinker.Tab)cell.agents[0];
				neighborSequenceList.setAsHeader(cellHeader, phase.speciesMaster.atomList.size());
				listIterator.setBasis(neighborSequenceList);
			}
		} else {
			listIterator.setBasis(a);
		}
	}
    
	/**
	 * Resets iterator so that it will loop up the list of atoms beginning
	 * from the first one.
	 */
	public Atom reset() {return listIterator.reset();}
    
	/**
	 * Sets to state in which hasNext is false.
	 */
	public void unset() {listIterator.unset();}
        
	public boolean hasNext() {return listIterator.hasNext();}
    
	public boolean contains(Atom atom) {return listIterator.contains(atom);}
    
	/**
	 * Resets for iteration according to the given directive.  If the directive does
	 * not specify an atom, this is the same as the reset() method, except that the
	 * direction of iteration is as given by the directive.  If an atom is specified,
	 * iteration begins with it and proceeds up or down list from there.
	 */
	public Atom reset(IteratorDirective id) {return listIterator.reset(id);}
    
	/**
	 * Returns the next atom in the iteration sequence.  Assumes that hasNext is
	 * true; calling when hasNext is false can lead to unpredictable results, and
	 * may or may not cause an error or exception.
	 */
	public Atom next() {
//		Atom myNext = listIterator.next(); System.out.println(myNext); return myNext;
		return listIterator.next();
	}
    
	public void allAtoms(AtomAction act) {
		listIterator.allAtoms(act);
	}
    
	public Atom getBasis() {
		return basis.atom;
	}
    
	public int size() {
		return listIterator.size();
	}
    
	/**
	 * Method to test SequentialIterator
	 */
	public static void main(String args[]) {
		Simulation sim = new Simulation(new Space2D());

		IteratorFactoryCell iteratorFactory = new IteratorFactoryCellLeaf(sim);
		sim.setIteratorFactory(iteratorFactory);
		int nAtoms = 6;       
		SpeciesSpheresMono speciesSpheres = new SpeciesSpheresMono();
		speciesSpheres.setNMolecules(nAtoms);
		Potential potential = new P2HardSphere();
		Phase phase = new Phase();
		IntegratorHard integrator = new IntegratorHard();
		integrator.setTimeStep(0.01);
		Controller controller = new Controller();
		Simulation.instance.elementCoordinator.go();
		System.out.println("Starting MD");
		integrator.setMaxSteps(100);
		integrator.initialize();
		integrator.run();
		System.out.println("Done");
		
		
		BravaisLattice lattice = ((IteratorFactoryCell)sim.getIteratorFactory()).getLattice(phase);
		AtomIterator iterator = iteratorFactory.makeGroupIteratorSequential();
		iterator.setBasis(phase.getAgent(speciesSpheres));
		Atom first = null;
		Atom middle = null;
		Atom last = null;
		iterator.reset();
		int k = 0;
		while(iterator.hasNext()) {
			Atom a = iterator.next();
			if(k == 0) first = a;
			else if(k == nAtoms/2) middle = a;
			else if(k == nAtoms-1) last = a;
			k++;
		}
	    
		IteratorDirective.testSuite(iterator, first, middle, last);
	    
	}//end of SequentialIterator.main
    
}//end of SequentialIterator
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Iterates neighbors from among all leaf atoms.
 */
public static final class IntragroupNbrIterator extends AtomIterator {
	
	private final IteratorDirective upSkip = new IteratorDirective(IteratorDirective.UP);
	private final IteratorDirective dnSkip = new IteratorDirective(IteratorDirective.DOWN);
    
	public IntragroupNbrIterator(IteratorFactoryCell factory) {
		iteratorFactory = factory;
		upSkip.setSkipFirst(true);
		dnSkip.setSkipFirst(true);
	}
    
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}
    
    /**
     * Atom basis is ignored
     * @see etomica.AtomIterator#all(Atom, IteratorDirective, AtomActive)
     */
	public void all(Atom dummy, IteratorDirective id, final AtomActive action) {
		if(action == null) return;
		Atom atom = id.atom1();
		if(atom == null) return;
		AbstractCell referenceCell = ((CellSequencer)atom.seq).cell();//cell in which reference atom resides
		boolean upListNow = id.direction().doUp();
		boolean doGoDown = id.direction().doDown();
		if(referenceCell != null) {
			AtomLinker tab = referenceCell.neighborManager().tab;
			if(upListNow) {
				//loop over rest of atoms in cell
				AtomIteratorList.all(((NeighborSequencer)atom.seq).nbrLink, null, upSkip, action);
				//loops over atoms in uplist cells
				for(AtomLinker e=tab.next; true; e=e.next) {//loop over cells
					if(e.atom == null && ((AtomLinker.Tab)e).isHeader()) break;
					AtomLinker.Tab tabs = (AtomLinker.Tab)e.atom.agents[0];
					AtomLinker next = tabs.next;
					while(next.atom != null) {//loop up atoms in cell
						action.actionPerformed(next.atom);
						next = next.next;
					}//end while
				}
			}//end if(upListNow)
			if(doGoDown) {
				//loop over rest of atoms in cell
				AtomIteratorList.all(((NeighborSequencer)atom.seq).nbrLink, null, dnSkip, action);
				//loops over atoms in dnlist cells
				for(AtomLinker e=tab.previous; true; e=e.previous) {//loop over cells
					if(e.atom == null && ((AtomLinker.Tab)e).isHeader()) break;
					AtomLinker.Tab tabs = (AtomLinker.Tab)e.atom.agents[0];
					AtomLinker next = tabs.nextTab.previous;//get the last atom on the cell marked by tabs[tabIndex]
					while(next.atom != null) {//loop down atoms in cell
						action.actionPerformed(next.atom);
						next = next.previous;
					}//end while
				}				
			}//end if(doGoDown)			
		} else { //no cell-based iteration
			AtomIteratorList.all(atom.seq, id, action);
		}
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
		return atom.node.isLeaf();
	}
    
	/**
	 * Sets to state in which hasNext is false.
	 */
	public void unset() {listIterator.unset();}
    
	/**
	 * Does reset if atom in iterator directive is child of the current basis.  
	 * Sets hasNext false if given atom does is not child of basis.  Throws
	 * an IllegalArgumentException if directive does not specify an atom.
	 */
	public Atom reset(IteratorDirective id) {
		direction = id.direction();
		return doReset(id.atom1());
	}
    
	//we assume that the only Tab links in the list are those demarking
	//the beginning of each cell's sequence; thus we reset the list iterator
	//using null as the terminator
    
	public Atom reset(Atom atom) {
		return doReset(atom);
	}
    
	private Atom doReset(Atom atom) {
		referenceAtom = atom;
		if(atom == null) 
			throw new IllegalArgumentException("Cannot reset IteratorFactoryCell.IntragroupNbrIterator without referencing an atom");

		referenceCell = ((NeighborSequencer)atom.seq).cell();
		cellIterator.setBasis(referenceCell);
		return doReset();
	}
    
	private Atom doReset() {
		listIterator.setSkipFirstAtom(true);
		upListNow = direction.doUp();
		doGoDown = direction.doDown();
		if(iterateCells) {
			listIterator.unset();
			if(upListNow) {
				cellIterator.reset(IteratorDirective.UP);//set cell iterator to return first up-neighbor of reference cell
				listIterator.reset(((NeighborSequencer)referenceAtom.seq).nbrLink, null, IteratorDirective.UP);
			}
			if(!listIterator.hasNext()) advanceCell();
		} else {//no cell iteration
			listIterator.reset(referenceAtom.seq, direction);
		}
		return listIterator.peek();
	}
	// Moves to next cell that has an iterate
	private void advanceCell() {
		do {
			if(cellIterator.hasNext()) {
				listIterator.setSkipFirstAtom(false);
				Atom cell = cellIterator.next();
				AtomLinker.Tab tabs = (AtomLinker.Tab)cell.agents[0];
				if(upListNow) {
					listIterator.reset(tabs, null, IteratorDirective.UP);
				} else {
					listIterator.reset(tabs.nextTab, null, IteratorDirective.DOWN);
				}
			} else if(doGoDown) {//no more cells that way; see if should now reset to look at down-cells
				listIterator.setSkipFirstAtom(true);
				cellIterator.reset(IteratorDirective.DOWN);//set cell iterator to return first down neighbor of reference cell
				listIterator.reset(((NeighborSequencer)referenceAtom.seq).nbrLink, null, IteratorDirective.DOWN);
				upListNow = false;
				doGoDown = false;
			} else {//no more cells at all
				break;
			}
		} while(!listIterator.hasNext());
	}

	/**
	 * Performs given action for each child atom of basis.
	 */
	public void allAtoms(AtomAction act) {
		doReset();
		if(iterateCells) {
			while(listIterator.hasNext()) {
				listIterator.allAtoms(act);
				advanceCell();
			}
		} else {
			listIterator.allAtoms(act);
		}
	}
            
	public Atom next() {
		Atom atom = listIterator.next();
//		System.out.println("nbr"+atom);
		if(!listIterator.hasNext() && iterateCells) advanceCell();
		return atom;
	}
    
	/**
	 * Throws RuntimeException because this is a neighbor iterator, and must
	 * be reset with reference to an atom.
	 */
	public Atom reset() {
		throw new RuntimeException("Cannot reset IteratorFactoryCell.IntragroupNbrIterator without referencing an atom");
	}
    
	/**
	 * Not implemented, as basis is not relevant to neighbor iteration for the
	 * this leaf iterator.
	 */
	public void setBasis(Atom atom) {
		setBasis((AtomTreeNodeGroup)atom.node);
	}
    
	/**
	 * Not implemented, as basis is not relevant to neighbor iteration for the
	 * this leaf iterator.
	 */
	public void setBasis(AtomTreeNodeGroup node) {
		//throw new RuntimeException("Unexpected attempt to set basis in neighbor iterator");
	}
    
	/**
	 * Returns null.
	 */
	public Atom getBasis() {return null;}
    
	/**
	 * The number of atoms in the current basis.  This will differ from
	 * the number of iterates given by the iterator, because the iterator
	 * will return only those atoms neighboring a reference atom, which
	 * in general is not all the atoms in the basis.
	 */
	public int size() {
	   throw new etomica.exception.MethodNotImplementedException();
	}   

	private Atom next;
	private Atom referenceAtom;
	private boolean upListNow, doGoDown;
	private IteratorDirective.Direction direction = IteratorDirective.BOTH;
	private AbstractCell referenceCell;
	private boolean iterateCells = true;
	private int tabIndex;
	private BravaisLattice lattice;
	private final SiteIteratorNeighbor cellIterator = new SiteIteratorNeighbor();
	private final AtomIteratorList listIterator = new AtomIteratorList();
	private final IteratorFactoryCell iteratorFactory;

	/**
	 * Method to test IntragroupNbrIterator
	 */
	public static void main(String args[]) {
		Default.ATOM_SIZE = 1.0;
		Simulation sim = new Simulation(new Space2D());

		IteratorFactoryCell iteratorFactory = new IteratorFactoryCell(sim);
		sim.setIteratorFactory(iteratorFactory);
		int nAtoms = 5;       
		SpeciesSpheresMono speciesSpheres = new SpeciesSpheresMono();
		speciesSpheres.setNMolecules(nAtoms);
		Potential potential = new P2HardSphere();
		Phase phase = new Phase();
		IntegratorHard integrator = new IntegratorHard();
		integrator.setTimeStep(0.01);
		Controller controller = new Controller();
		Simulation.instance.elementCoordinator.go();
		System.out.println("Starting MD");
		integrator.setMaxSteps(100);
		integrator.initialize();
		integrator.run();
		System.out.println("Done");
		
		
		BravaisLattice lattice = ((IteratorFactoryCell)sim.getIteratorFactory()).getLattice(phase);
		AtomIterator sequenceIterator = iteratorFactory.makeGroupIteratorSequential();
		AtomIterator iterator = iteratorFactory.makeIntragroupNbrIterator();
		sequenceIterator.setBasis(phase.getAgent(speciesSpheres));
		iterator.setBasis(phase.getAgent(speciesSpheres));
		Atom first = null;
		Atom middle = null;
		Atom last = null;
		sequenceIterator.reset();
		int k = 0;
		while(sequenceIterator.hasNext()) {
			Atom a = sequenceIterator.next();
			if(k == 0) first = a;
			else if(k == nAtoms/2) middle = a;
			else if(k == nAtoms-1) last = a;
			k++;
		}
		System.out.println(" first: "+first.signature());
		System.out.println("middle: "+middle.signature());
		System.out.println("  last: "+last.signature());
	    
		IteratorDirective.testSuite(iterator, first, middle, last);
	    
	}//end of IntragroupNbrIterator.main

}//end of IntragroupNbrIterator class


/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Sequencer used for atoms being neighbor listed.
 */

public static final class NeighborSequencer extends AtomSequencer implements CellSequencer {
    
	public AbstractCell cell;         //cell currently occupied by this atom
	public BravaisLattice lattice;    //cell lattice in the phase occupied by this atom
	public final AtomLinker nbrLink;  //linker used to arrange atom in sequence according to cells
    
	public NeighborSequencer(Atom a) {
		super(a);
		nbrLink = new AtomLinker(a);
	}
    
	//CellSequencer interface method
	public void latticeChangeNotify() {
		this.assignCell();
		if(atom.node.isLeaf()) return;
		else throw new RuntimeException("Unexpected occurrence of NeighborSequencer in a non-leaf atom");
	}

	//CellSequencer interface method
	public AbstractCell cell() {return cell;}

	public void remove() {
		super.remove();
		nbrLink.remove();
	}
        
	public void addBefore(AtomLinker newNext) {
		//newNext will never be one of the cell tabs
		super.addBefore(newNext);
		assignCell();
	}
	/**
	 * Reshuffles position of "neighbor" links without altering the regular links.
	 */
	public void moveBefore(AtomLinker newNext) {
		nbrLink.moveBefore(newNext);
	}


	/**
	 * Returns true if this atom preceeds the given atom in the atom sequence.
	 * Returns false if the given atom is this atom, or (of course) if the
	 * given atom instead preceeds this one.
	 * 
	 * Not implemented
	 */
	public boolean preceeds(Atom a) {
		return atom.node.index() < a.node.index();
//		throw new etomica.exception.MethodNotImplementedException();
	}
    
	/**
	 * Method called when a translate method of coordinate is invoked.
	 */
	public void moveNotify() {
		assignCell();
	}
    
	/**
	 * Method called when the parent of the atom is changed.
	 * By the time this method is called, the atom has been placed
	 * in the childList of the given parent (if it is not null).
	 */
	public void setParentNotify(AtomTreeNodeGroup newParent) {
		if(newParent == null || newParent instanceof AtomReservoir.ReservoirAtomTreeNode) {
			lattice = null;
		} else {
			Simulation sim = newParent.parentSimulation();
			if(sim == null) lattice = null;
			else lattice = ((IteratorFactoryCell)sim.iteratorFactory).getLattice(newParent.parentPhase());
		}
//		if(lattice != null) {
//			Atom testCell = lattice.siteList().getFirst();
//			if(testCell.agents == null) setupListTabs(lattice);
//		}
		assignCell();
	}

//Determines appropriate cell and assigns it
	public void assignCell() {
		if(lattice == null) {
			assignCell(null);
		} else {
			int[] latticeIndex = lattice.getPrimitive().latticeIndex(atom.coord.position(), lattice.getDimensions());
			AbstractCell newCell = (AbstractCell)lattice.site(latticeIndex);
			if(newCell != cell) {assignCell(newCell);}
		}
	}
//Assigns atom to given cell
	public void assignCell(AbstractCell newCell) {
		cell = newCell;
		if(newCell == null) {
			nbrLink.remove();
		} else {
			this.moveBefore(((AtomLinker.Tab)newCell.agents[0]).nextTab);
		}
	}//end of assignCell
    
	public static final AtomSequencer.Factory FACTORY = new AtomSequencer.Factory() {
		public AtomSequencer makeSequencer(Atom atom) {return new NeighborSequencer(atom);}
		public Class sequencerClass() {return NeighborSequencer.class;}
	};
}//end of NeighborSequencer

/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Iterates through the cells of the given lattice, and adds a tab to the 
 * given list for each cell, and extends the tablist in each cell to reference
 * its new tab.
 */
private static void setupListTabs(BravaisLattice lattice) {
	AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
	iterator.reset();
	AtomLinker.Tab header = null;
	while(iterator.hasNext()) {//loop over cells
		Site site = (Site)iterator.next();
		if(site.agents == null) site.agents = new Object[1];
		AtomLinker.Tab newTab = AtomLinker.newTab();
		if(header == null) header = newTab;
		else newTab.addBefore(header);
		site.agents[0] = newTab;
	}
}//end of setupListTabs

	/**
	 * Demonstrates how this class is implemented.
	 */
	public static void main(String[] args) {
        
  //      SequentialIterator.main(args); 
  //      IntragroupNbrIterator.main(args); 
  //      IntergroupNbrIterator.main(args); 
		boolean cellListing = true;
		boolean doDisplay = true;
		boolean md = true;
		int nAtoms = 90;
		boolean fixedLengthRun = false;
		boolean mixture = true;
        
		System.out.println(cellListing ? "cellListing" : "no cellListing");
		System.out.println(md ? "md" : "mc");
		System.out.println("number of atoms in simulation: "+2*nAtoms);
        
		Default.ATOM_SIZE = 150./(double)nAtoms;
//		Default.ATOM_SIZE = 3;
		etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space2D());
		Simulation.instance = sim;

		if(cellListing) sim.setIteratorFactory(new IteratorFactoryCellLeaf(sim));
        
		Integrator integrator = null;
		if(md) {
			integrator = new IntegratorHard();
			((IntegratorHard)integrator).setTimeStep(0.01);
		} else {
			integrator = new IntegratorMC();
			MCMoveAtom mcMoveAtom = new MCMoveAtom((IntegratorMC)integrator);
			//MCMoveVolume mcMoveVolume = new MCMoveVolume((IntegratorMC)integrator);
			//mcMoveVolume.setPressure(3.0);
		}
        
		Controller controller = new Controller();
		etomica.graphics.DisplayPhase displayPhase = null;
		if(doDisplay) displayPhase = new etomica.graphics.DisplayPhase();

		if(!mixture) {
			SpeciesSpheresMono speciesSpheres = new SpeciesSpheresMono();
			speciesSpheres.setNMolecules(2*nAtoms);
			if(doDisplay) etomica.graphics.ColorSchemeByType.setColor(speciesSpheres, java.awt.Color.green);
	    	    
			Potential2 potential = new P2HardSphere();
			potential.setSpecies(speciesSpheres);
//			potential.setIterator(new ApiGeneral(sim.space,sim.iteratorFactory.makeGroupIteratorSequential(),
//												  sim.iteratorFactory.makeIntragroupNbrIterator()));
		} else {
			SpeciesSpheresMono speciesSpheres1 = new SpeciesSpheresMono();
			SpeciesSpheresMono speciesSpheres2 = new SpeciesSpheresMono();
			speciesSpheres1.setNMolecules(nAtoms);
			speciesSpheres2.setNMolecules(nAtoms);
			if(doDisplay) {
				etomica.graphics.ColorSchemeByType.setColor(speciesSpheres1, java.awt.Color.green);
				etomica.graphics.ColorSchemeByType.setColor(speciesSpheres2, java.awt.Color.black);
			}
	    	    
			Potential2 potential11 = new P2HardSphere();
		 //   Potential2 potential12 = new P2HardSphere();
//			Potential2 potential12 = new P2SquareWell();
//			Potential2 potential22 = new P2HardSphere();
			potential11.setSpecies(speciesSpheres1, speciesSpheres1);
//			potential12.setSpecies(speciesSpheres2, speciesSpheres1);
//			potential22.setSpecies(speciesSpheres2, speciesSpheres2);
		}
        
		Phase phase = new Phase();//must come after species!
//		etomica.graphics.DeviceTrioControllerButton buttons = new etomica.graphics.DeviceTrioControllerButton(controller);
       
		etomica.graphics.LatticeRenderer.ColorSchemeNeighbor colorSchemeNbr = null; 
		etomica.graphics.LatticeRenderer.ColorSchemeCell colorSchemeCell = null;
		if(doDisplay) {
			if(cellListing) {
//				colorSchemeNbr = new etomica.graphics.LatticeRenderer.ColorSchemeNeighbor(sim);
				colorSchemeCell = new etomica.graphics.LatticeRenderer.ColorSchemeCell();
			}
		}
    	    
		integrator.setDoSleep(true);
		integrator.setSleepPeriod(2);
	    
	    
		//this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		
		if(cellListing) {

			BravaisLattice lattice = ((IteratorFactoryCell)sim.getIteratorFactory()).getLattice(phase);
		//    ((IteratorFactoryCell)sim.iteratorFactory).setNeighborRange(4.0);
    	
			if(doDisplay) {
				//draw lattice cells on display
		//      etomica.graphics.LatticeRenderer latticeRenderer = 
		//              new etomica.graphics.LatticeRenderer(lattice);
		//      displayPhase.addDrawable(latticeRenderer);
    	  
				//color atoms
//				colorSchemeNbr.setAtom(phase.speciesMaster.atomList.getRandom());
//				displayPhase.setColorScheme(colorSchemeNbr);
		//        colorSchemeCell.setLattice(lattice);
		//	    displayPhase.setColorScheme(colorSchemeCell);
			}
		}
        
		if(doDisplay) etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
        
		if(fixedLengthRun) {
			if(md) integrator.setMaxSteps(1000);
			else   integrator.setMaxSteps(100*nAtoms*2);
		}
//		  integrator.initialize();
		System.out.println("Starting");
//		  etomica.benchmark.Stopwatch timer = new etomica.benchmark.Stopwatch().start();
//		  integrator.run();
//		  timer.stop();
//		  System.out.println("Elapsed time: "+0.001*timer.getElapsedTime());
        
	 //   controller.start();
	}//end of main
    
    
   
}//end of IteratorFactoryCell