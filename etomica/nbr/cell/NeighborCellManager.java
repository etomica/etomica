/*
 * History
 * Created on Nov 21, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.AtomIteratorList;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorMolecule;
import etomica.AtomLinker;
import etomica.AtomList;
import etomica.Default;
import etomica.Phase;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.Space;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeEvent;
import etomica.lattice.LatticeListener;
import etomica.lattice.Primitive;
import etomica.lattice.PrimitiveCubic;
import etomica.lattice.PrimitiveOrthorhombic;
import etomica.lattice.Site;
import etomica.math.geometry.Polyhedron;
import etomica.nbr.cell.IteratorFactoryCell.CellSequencer;
import etomica.utility.Arrays;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class NeighborCellManager {

    private final BravaisLattice lattice;
    private final int[] dimensions;
    private final Space space;
    private final Primitive primitive;
    private final Phase phase;
    protected double neighborRange;
    private int listCount;
    
    public NeighborCellManager(Phase phase, int nCells) {
        this(phase, nCells, new PrimitiveCubic(phase.space()));
    }
    /**
     * 
     */
    public NeighborCellManager(Phase phase, int nCells, Primitive primitive) {
        this.phase = phase;
        space = phase.space();
        this.primitive = primitive;
        neighborRange = Default.ATOM_SIZE;
        dimensions = new int[space.D()];
        for(int i=0; i<space.D(); i++) dimensions[i] = nCells;
        lattice = makeCellLattice(phase);
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
    private BravaisLattice makeCellLattice(final Phase phase) {
        //make the unit cell factory and set it to produce cells of the appropriate size
        final Primitive primitiveCopy = primitive.copy();//each new lattice works with its own primitive
        Space.Vector primitiveSize = space.makeVector();
        primitiveSize.E(phase.boundary().dimensions());
        primitiveSize.DE(Space.makeVector(dimensions));
        primitiveCopy.setSize(primitiveSize.toArray());
        //construct the lattice
        final BravaisLattice lattice = 
            BravaisLattice.makeLattice(space, NeighborCell.makeFactory(space,primitiveCopy),
                                        dimensions, primitiveCopy);
        lattice.shiftFirstToOrigin();
        primitiveCopy.setLattice(null);
        
        //set up the neighbor lists for each cell in the lattice
        //TODO merge with etomica.nbr.NeighborCriterion, or just use AtomFilter
        etomica.lattice.NeighborManager.Criterion neighborCriterion = new etomica.lattice.NeighborManager.Criterion() {
            final Space.Vector dr = space.makeVector();
            final Space.CoordinatePair cPair;
            //initializer
            {   cPair = space.makeCoordinatePair();
                cPair.setBoundary(phase.boundary());
             }
            public boolean areNeighbors(Site s1, Site s2) {
                cPair.reset(s1.coord, s2.coord);
                return ((AtomTypeCell)s1.type).cell.r2NearestVertex(cPair.dr()) < neighborRange*neighborRange;
            }
        };
        lattice.setupNeighbors(neighborCriterion);
                
        //add listener to notify all sequencers of any lattice events (resizing of lattice, for example)
        lattice.eventManager.addListener(new LatticeListener() {
            public void actionPerformed(SimulationEvent evt) {
                actionPerformed((LatticeEvent)evt);
            }
            public void actionPerformed(LatticeEvent evt) {
                if(evt.type() == LatticeEvent.REBUILD || evt.type() == LatticeEvent.ALL_SITE) {
                    AtomIteratorMolecule iterator = new AtomIteratorMolecule(phase);
                    iterator.reset();
                    while(iterator.hasNext()) {
                        ((CellSequencer)iterator.nextAtom().seq).latticeChangeNotify();
                    }
                }//end if
            }
        });
        
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
        return lattice;
    }//end of makeCellLattice method
            
    /**
     * Sets the maximum range of interaction for which the cells must keep neighbors.
     * Updates the neighbor lists of all previous and future cell lattices consistent
     * with the new value.
     */
    public void setNeighborRange(double r) {
        neighborRange = r;
        final double rangeSquared = r*r;
        etomica.lattice.NeighborManager.Criterion neighborCriterion = new etomica.lattice.NeighborManager.Criterion() {
            public boolean areNeighbors(Site s1, Site s2) {
                return ((Polyhedron)s1).r2NearestVertex((Polyhedron)s2, phase.boundary()) < rangeSquared;
            }
        };
        lattice.setupNeighbors(neighborCriterion);
    }
    /**
     * Accessor method for the maximum range of interaction for which 
     * cells must keep neighbors.
     */
    public double getNeighborRange() {return neighborRange;}

    /**
     * @return the number of atom lists held by each cell
     */
    public int getListCount() {
        return listCount;
    }
    
    public BravaisLattice getCellLattice() {
        return lattice;
    }
    
    /**
     * Adds an atom list to each cell of the lattice.
     */
    public void addList() {
        AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
        iterator.reset();
        listCount++;
        while(iterator.hasNext()) {//loop over cells
            Site site = (Site)iterator.nextAtom();
            if(site.agents == null || site.agents[0] == null) {
                site.agents = new Object[1];
                site.agents[0] = new AtomList[1];
                ((AtomList[])site.agents[0])[0] = new AtomList();
                return;
            }
            Arrays.addObject((Object[])site.agents[0],new AtomList());
        }
    }//end of setupListTabs
    


}
