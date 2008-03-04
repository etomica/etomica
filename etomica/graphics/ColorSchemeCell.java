package etomica.graphics;

import java.awt.Color;
import java.util.HashMap;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.lattice.FiniteLattice;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.NeighborCellManager;
import etomica.box.BoxAgentManager;

public class ColorSchemeCell extends ColorScheme {
    
    public ColorSchemeCell(PotentialMasterNbr potentialMaster, IRandom random, IBox box) {
        BoxAgentManager cellAgentManager = potentialMaster.getCellAgentManager();
        cellManager = (NeighborCellManager)cellAgentManager.getAgent(box);
        this.random = random;
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
