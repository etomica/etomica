/*
 * History
 * Created on Nov 21, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomIteratorList;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorMolecule;
import etomica.AtomIteratorPhaseDependent;
import etomica.Default;
import etomica.Integrator;
import etomica.Phase;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.Space;
import etomica.SpeciesAgent;
import etomica.lattice.CellLattice;
import etomica.lattice.PrimitiveCubic;
import etomica.lattice.PrimitiveOrthorhombic;

/**
 * Class that defines and manages construction and use of lattice of cells 
 * for cell-based neighbor listing.
 */

public class NeighborCellManager implements Integrator.IntervalListener {

    private final CellLattice lattice;
    private final int[] dimensions;
    private final Space space;
    private final Phase phase;
    protected double neighborRange;
    private int listCount;
    private final AtomIteratorPhaseDependent atomIterator;
    private int iieCount;
    private int updateInterval;
    private int priority;
    
    /**
     * 
     */
    public NeighborCellManager(Phase phase, int nCells) {
        this.phase = phase;
        space = phase.space();
        neighborRange = Default.ATOM_SIZE;
        dimensions = new int[space.D()];
        for(int i=0; i<space.D(); i++) dimensions[i] = nCells;
        lattice = makeCellLattice(phase);
        atomIterator = new AtomIteratorMolecule(phase);
        setPriority(150);
        setUpdateInterval(1);
        addList(phase.speciesMaster.node.childAtomCount());//add occupant lists to cells for each species already present in phase
    }
    /**
     * Constructs the cell lattice used to organize all cell-listed atoms
     * in the given phase.  Each new lattice is set up with its own instance
     * of the primitive, formed by copying the instance associated with this factory.  
     * Note that the phase does not contain any reference
     * to the lattice.  Its association with the phase is made through the 
     * deployedLattices array kept by this iterator factory class, and by 
     * the reference in each neighbor sequencer to the cell containing its atom.
     */
    private CellLattice makeCellLattice(final Phase phase) {
        //construct the lattice
        final CellLattice lattice = new CellLattice(space.D, NeighborCell.FACTORY);
        Space.Vector cellSize = (Space.Vector)phase.boundary().dimensions().clone();
        cellSize.DE(Space.makeVector(dimensions));
        lattice.setCellSize(cellSize.toArray());
                
        //add listener to phase to update the size and placement of the lattice
        //cells if the phase undergoes an inflation of its boundary.
        //An inflation event should not, however, cause the molecules to be reassigned
        //to their lattice cells, since the atom positions and the cells scale proportionately
        //Atoms in molecules may be reassigned, if they are the focus of neighbor listing (which is not usually the case)
        phase.boundaryEventManager.addListener(new PhaseListener() {
            final AtomIteratorListSimple cellIterator = new AtomIteratorListSimple();
            double[] newSize;
            public void actionPerformed(SimulationEvent evt) {
                actionPerformed((PhaseEvent)evt);
            }
            public void actionPerformed(PhaseEvent evt) {
                if(!evt.type().equals(PhaseEvent.BOUNDARY_INFLATE)) return;
                    //we expect that primitive.lattice() == null, so change of size doesn't cause replacement of atoms in cells
                if(evt.isotropic) {
                    primitiveCopy.scaleSize(evt.isoScale);
                    cellIterator.setList(lattice.siteList());
                    cellIterator.reset();
                    while(cellIterator.hasNext()) {
                        cellIterator.nextAtom().coord.inflate(evt.isoScale);
                    }
                } else {//anisotropic inflation
                    if(primitiveCopy instanceof PrimitiveCubic) throw new RuntimeException("Cannot handle anisotropic inflate with cubic primitive used in IteratorFactoryCell");
                    newSize = ((PrimitiveOrthorhombic)primitiveCopy).getSize();
                    for(int i=0; i<newSize.length; i++) newSize[i] *= evt.anisoScale.x(i);
                    primitiveCopy.setSize(newSize);
                    cellIterator.setList(lattice.siteList());
                    cellIterator.reset();
                    while(cellIterator.hasNext()) {
                        cellIterator.nextAtom().coord.inflate(evt.anisoScale);
                    }
                }
            }//end actionPerformed
        });
        
        //listener to phase to detect addition of new SpeciesAgent
        //or new atom
        phase.speciesMaster.addListener(new PhaseListener() {
            public void actionPerformed(SimulationEvent evt) {
                actionPerformed((PhaseEvent)evt);
            }
           public void actionPerformed(PhaseEvent evt) {
                if(evt.type() == PhaseEvent.ATOM_ADDED) {
                    Atom atom = evt.atom();
                    //new species agent requires another list in each cell
                    if(atom instanceof SpeciesAgent) {
                        addList();
                   //otherwise new atom placed in cell if at molecule level
                    } else if(atom.node.depth() == 2) {
                        assignCell(atom);
                    }
                }
            }
        });
        
        return lattice;
    }//end of makeCellLattice method
            
    /**
     * Sets the maximum range of interaction for which the cells must keep neighbors.
     * Updates the neighbor lists of all previous and future cell lattices consistent
     * with the new value.
     */
    public void setNeighborRange(double r) {
        neighborRange = r;
        lattice.setupNeighbors(makeNeighborCriterion());
    }
    /**
     * Accessor method for the maximum range of interaction for which 
     * cells must keep neighbors.
     */
    public double getNeighborRange() {return neighborRange;}

    /**
     * @return the number of atom lists held by each cell.
     */
    public int getListCount() {
        return listCount;
    }
    
    /**
     * @return the lattice of cells.
     */
    public CellLattice getCellLattice() {
        return lattice;
    }
    
    /**
     * Assigns cells to all molecules in the phase.
     */
    public void assignCellAll() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            assignCell(atomIterator.nextAtom());
        }
    }
    
    /**
     * Assigns the cell for the given atom.
     * @param atom
     */
    public void assignCell(Atom atom) {
        AtomSequencerCell seq = (AtomSequencerCell)atom.seq;
        NeighborCell newCell = (NeighborCell)lattice.site(atom.coord.position());
        if(newCell != seq.cell) {assignCell(seq, newCell, atom.node.parentSpeciesAgent().node.index());}
    }
    
    /**
     * Assigns atom sequencer to given cell in the list of the given index.
     */
    public void assignCell(AtomSequencerCell seq, NeighborCell newCell, int listIndex) {
        if(seq.cell != null) seq.cell.occupants()[listIndex].remove(seq.nbrLink);
        seq.cell = newCell;
//        seq.nbrLink.remove();
        if(newCell != null) {
            newCell.occupants()[listIndex].add(seq.nbrLink);
        }
    }//end of assignCell
    
    private void addList(int n) {
        for(int i=0; i<n; i++) addList();
    }
    
    /**
     * Adds an AtomList to each cell of the lattice.  This is performed
     * when a new species is added to the simulation.  Each list associated
     * with a cell holds the molecules of a given species that are in that cell.
     */
    public void addList() {
        listCount++;
        Object[] sites = lattice.sites();
        for(int i=sites.length-1; i>=0; i--) {
            ((NeighborCell)sites[i]).addOccupantList();
        }
    }

    public void intervalAction(Integrator.IntervalEvent event) {
        if (event.type() == Integrator.IntervalEvent.INITIALIZE) {
            assignCellAll();
        }
        else if (event.type() == Integrator.IntervalEvent.INTERVAL) {
            if (--iieCount == 0) {
                assignCellAll();
                iieCount = updateInterval;
            }
        }
    }
    
    public void setUpdateInterval(int updateInterval) {
        this.updateInterval = updateInterval;
    }

    public int getUpdateInterval() {
        return updateInterval;
    }
    
    public int getPriority() {
        return priority;
    }
    
    public void setPriority(int priority) {
        this.priority = priority;
    }
}
