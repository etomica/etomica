/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.ApiBuilder;
import etomica.ApiInnerFixed;
import etomica.Atom;
import etomica.AtomIteratorListSimple;
import etomica.AtomsetIteratorPhaseDependent;
import etomica.AtomsetIteratorTargetable;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.lattice.CellLattice;
import etomica.lattice.SimpleLattice;
import etomica.lattice.SiteIterator;

/**
 * Returns iterates formed from from the childList of two species.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the directive atom is descended.  
 */
public class ApiIntergroupCell implements
        AtomsetIteratorPhaseDependent, AtomsetIteratorTargetable {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
     * @param D the dimension of the space of the simulation (used to construct cell iterator)
	 */
	public ApiIntergroupCell(int D, Species[] species) {
        listIterator = ApiBuilder.makeInterlistIterator();
        cellNbrIterator = new CellLattice.NeighborIterator(D);
        aiInner = ((AtomIteratorListSimple)listIterator.getInnerIterator());
        aiOuter = ((AtomIteratorListSimple)listIterator.getOuterIterator());
        this.species = (Species[])species.clone();
        index0 = species[0].getIndex();
        index1 = species[1].getIndex();

	}

	public void setTarget(Atom[] targetAtoms) {
	    if(targetAtoms[0] == null) {
            cellIterator = allCellIterator;
        }
	}

	public void setPhase(Phase phase) {
		cellNbrIterator.setLattice(phase.getLattice());
        speciesAgent0 = phase.getAgent(species[0]);
        speciesAgent1 = phase.getAgent(species[1]);
        unset();
	}
    
    public void allAtoms(AtomsetAction action) {
        cellIterator.reset();
        while(cellIterator.hasNext()) {
            NeighborCell cell = (NeighborCell)cellIterator.next();
            aiOuter.setList(cell.occupants()[index0]);
            if(aiOuter.size() == 0) continue;//cell is empty of specie0 molecules
            //molecules in same cell
            aiInner.setList(cell.occupants()[index1]);
            if(aiInner.size() > 0) listIterator.allAtoms(action);
            //loop over neighbor cells
            cellNbrIterator.setSite(index);
            cellNbrIterator.reset();
            while(cellNbrIterator.hasNext()) {
                cell = (NeighborCell)cellIterator.next();
                aiInner.setList(cell.occupants()[index1]);
                if(aiInner.size() > 0) listIterator.allAtoms(action);
            }
        }
        
    }
	
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
	}
	
	/**
	 * Indicates whether the given atom pair will be among the iterates
	 * given by the iterator if reset in its present state.  True only
	 * if an iterated pair would match the atoms as ordered in the given
	 * array.
	 */
	public boolean contains(Atom[] atoms) {
        if(atoms==null || atoms[0]==null || atoms[1]==null || atoms[0]==atoms[1]) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

    public boolean hasNext() {
        // TODO Auto-generated method stub
        return false;
    }
    
    public Atom[] next() {
        // TODO Auto-generated method stub
        return null;
    }
    public Atom[] peek() {
        // TODO Auto-generated method stub
        return null;
    }
    public void reset() {
        // TODO Auto-generated method stub

    }
    public void unset() {
        // TODO Auto-generated method stub

    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    private final ApiInnerFixed listIterator;
    private final CellLattice.NeighborIterator cellNbrIterator;
    private SimpleLattice.Iterator allCellIterator;
    private SiteIterator cellIterator;
    private final Species[] species;
    private SpeciesAgent speciesAgent0, speciesAgent1;
    private final int index0, index1;
    private final AtomIteratorListSimple aiInner, aiOuter;
}
