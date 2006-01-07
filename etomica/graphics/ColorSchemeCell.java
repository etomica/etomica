package etomica.graphics;

import java.awt.Color;
import java.util.HashMap;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.FiniteLattice;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.phase.Phase;

/*
 * History
 * Created on Aug 23, 2005 by kofke
 */
public class ColorSchemeCell extends ColorSchemeCollective {
    
    private final HashMap hash = new HashMap();
    private final AtomIteratorListSimple allIterator = new AtomIteratorListSimple();
    
    public ColorSchemeCell(PotentialMasterCell potentialMasterCell) {
        super();
        this.potentialMasterCell = potentialMasterCell;
    }
    
    public void setLattice(FiniteLattice lattice) {
        Object[] sites = lattice.sites();
        for(int i=0; i<sites.length; i++) {
            hash.put(sites[i], ConstantsGraphic.randomColor());
        }
    }
    
    public void colorAllAtoms(Phase phase) {
        NeighborCellManager cellManager = potentialMasterCell.getNbrCellManager(phase);
        allIterator.setList(phase.getSpeciesMaster().atomList);
        allIterator.reset();
        while(allIterator.hasNext()) {
            Atom atom = allIterator.nextAtom();
            Object cell = cellManager.getCell(atom);
            atomColors[atom.getGlobalIndex()] = (Color)hash.get(cell);
        }
    }
    
    private final PotentialMasterCell potentialMasterCell;
}