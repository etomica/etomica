package etomica.graphics;

import java.awt.Color;
import java.util.HashMap;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.lattice.FiniteLattice;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

/*
 * History
 * Created on Aug 23, 2005 by kofke
 */
public class ColorSchemeCell extends ColorSchemeCollective {
    
    private final HashMap hash = new HashMap();
    private final AtomIteratorLeafAtoms allIterator = new AtomIteratorLeafAtoms();
    
    public ColorSchemeCell(Simulation sim) {
        super(sim);
        potentialMasterCell = (PotentialMasterCell)sim.potentialMaster;
    }
    
    public void setLattice(FiniteLattice lattice) {
        Object[] sites = lattice.sites();
        for(int i=0; i<sites.length; i++) {
            hash.put(sites[i], ConstantsGraphic.randomColor());
        }
    }
    
    public void colorAllAtoms(Phase phase) {
        NeighborCellManager cellManager = potentialMasterCell.getNbrCellManager(phase);
        allIterator.setPhase(phase);
        allIterator.reset();
        while(allIterator.hasNext()) {
            Atom atom = allIterator.nextAtom();
            Object cell = cellManager.getCell(atom);
            atomColors[atom.getGlobalIndex()] = (Color)hash.get(cell);
        }
    }
    
    private final PotentialMasterCell potentialMasterCell;
}