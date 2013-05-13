package etomica.nbr.cell;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.IAtomPositionDefinition;
import etomica.box.BoxAgentManager;
import etomica.nbr.site.PotentialMasterSite;
import etomica.space.ISpace;

/**
 * A PotentialMaster for use with a Simulation where cell-listing of atoms and
 * their neighbors is appropriate.
 * 
 * @author Andrew Schultz
 */
public class PotentialMasterCell extends PotentialMasterSite {

    /**
     * Creates PotentialMasterCell with default (1.0) range.  Range
     * should be set manually via setRange method.
     */
    public PotentialMasterCell(ISimulation sim, ISpace _space) {
        this(sim,1.0, _space);
    }
    
    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     * 
     * @param space the governing Space
     * @param range the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(ISimulation sim, double range, ISpace _space) {
        this(sim, range, (IAtomPositionDefinition)null, _space);
    }

    public PotentialMasterCell(ISimulation sim, double range,
            IAtomPositionDefinition positionDefinition, ISpace _space) {
        this(sim, range, new BoxAgentSourceCellManager(sim, positionDefinition, _space), _space);
    }
    
    public PotentialMasterCell(ISimulation sim, double range, 
    		BoxAgentSourceCellManager boxAgentSource, ISpace _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager(boxAgentSource), _space);
    }
    
    public PotentialMasterCell(ISimulation sim, double range, BoxAgentSourceCellManager boxAgentSource,
            BoxAgentManager agentManager, ISpace _space) {
        super(sim, boxAgentSource, agentManager, new Api1ACell(_space.D(),range,agentManager));
        setRange(range);
    }
    
    public double getRange() {
        return range;
    }
    
    public void setCellRange(int d) {
        super.setCellRange(d);

        BoxAgentManager.AgentIterator iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManager cellManager = (NeighborCellManager)iterator.next();
            cellManager.setCellRange(getCellRange());
        }
    }
    
    public void setRange(double d) {
        ((Api1ACell)neighborIterator).getNbrCellIterator().setNeighborDistance(d);
        ((BoxAgentSourceCellManager)boxAgentSource).setRange(d);
        range = d;

        BoxAgentManager.AgentIterator iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManager cellManager = (NeighborCellManager)iterator.next();
            cellManager.setPotentialRange(range);
        }
        
    }
    
    public NeighborCellManager getNbrCellManager(IBox box) {
        NeighborCellManager manager = (NeighborCellManager)boxAgentManager.getAgent(box);
        manager.setPotentialRange(range);
        int cr = getCellRange();
        if (cr == 0) throw new RuntimeException("need to set cell range first");
        manager.setCellRange(cr);
        return manager;
    }
    
    
    /**
     * Reassign atoms to cell lists for all boxes.
     */
    public void reset() {
        BoxAgentManager.AgentIterator iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManager neighborCellManager = (NeighborCellManager)iterator.next();
            neighborCellManager.assignCellAll();
        }
    }
    
    private static final long serialVersionUID = 1L;
    private double range;
}
