package etomica.graphics;

import java.util.HashMap;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.FiniteLattice;
import etomica.nbr.cell.AtomSequencerCell;
import etomica.phase.Phase;



/*
 * History
 * Created on Aug 23, 2005 by kofke
 */
public class ColorSchemeCell extends ColorSchemeCollective {
    
    private final HashMap hash = new HashMap();
    private final AtomIteratorListSimple allIterator = new AtomIteratorListSimple();
    
    public void setLattice(FiniteLattice lattice) {
        Object[] sites = lattice.sites();
        for(int i=0; i<sites.length; i++) {
            hash.put(sites[i], ConstantsGraphic.randomColor());
        }
    }
    
    public void colorAllAtoms(Phase phase) {
        allIterator.setList(phase.getSpeciesMaster().atomList);
        allIterator.reset();
        while(allIterator.hasNext()) {
            Atom atom = allIterator.nextAtom();
            Object cell = ((AtomSequencerCell)atom.seq).getCell();
            atom.allatomAgents[agentIndex] = hash.get(cell);
        }
    }
}