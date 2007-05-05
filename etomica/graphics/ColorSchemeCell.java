package etomica.graphics;

import java.awt.Color;
import java.util.HashMap;

import etomica.atom.IAtom;
import etomica.lattice.FiniteLattice;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.NeighborCellManager;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.simulation.Simulation;
import etomica.util.IRandom;

public class ColorSchemeCell extends ColorScheme {
    
    public ColorSchemeCell(Simulation sim, Phase phase) {
        PhaseAgentManager cellAgentManager = ((PotentialMasterNbr)sim.getPotentialMaster()).getCellAgentManager();
        cellManager = (NeighborCellManager)cellAgentManager.getAgent(phase);
        random = sim.getRandom();
    }
    
    public void setLattice(FiniteLattice lattice) {
        Object[] sites = lattice.sites();
        for(int i=0; i<sites.length; i++) {
            hash.put(sites[i], new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble()));
        }
    }
    
    public Color getAtomColor(IAtom a) {
        return (Color)hash.get(cellManager.getCell(a));
    }
    
    private static final long serialVersionUID = 1L;
    private final HashMap hash = new HashMap();
    private final NeighborCellManager cellManager;
    private final IRandom random;
}
