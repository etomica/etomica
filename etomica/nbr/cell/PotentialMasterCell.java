package etomica.nbr.cell;

import etomica.api.IBox;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomType;
import etomica.nbr.site.PotentialMasterSite;
import etomica.box.BoxAgentManager;
import etomica.potential.IPotential;
import etomica.potential.PotentialArray;
import etomica.simulation.ISimulation;
import etomica.space.Space;

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
    public PotentialMasterCell(ISimulation sim, Space _space) {
        this(sim,1.0, _space);
    }
    
    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     * 
     * @param space the governing Space
     * @param range the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(ISimulation sim, double range, Space _space) {
        this(sim, range, (AtomPositionDefinition)null, _space);
    }

    public PotentialMasterCell(ISimulation sim, double range,
            AtomPositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManager(positionDefinition, _space));
    }
    
    public PotentialMasterCell(ISimulation sim, double range, BoxAgentSourceCellManager boxAgentSource) {
        this(sim, range, boxAgentSource, new BoxAgentManager(boxAgentSource));
    }
    
    public PotentialMasterCell(ISimulation sim, double range, BoxAgentSourceCellManager boxAgentSource,
            BoxAgentManager agentManager) {
        super(sim, boxAgentSource, agentManager, new Api1ACell(sim.getSpace().D(),range,agentManager));
        setRange(range);
    }
    
    public double getRange() {
        return range;
    }
    
    /**
     * Returns the maximum range of any potential held by this potential master
     */
    public double getMaxPotentialRange() {
        return maxPotentialRange;
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
        manager.setCellRange(getCellRange());
        return manager;
    }
    
    
    /**
     * Recomputes the maximum potential range (which might change without this
     * class receiving notification) and readjust cell lists if the maximum
     * has changed.
     */
    public void reset() {
        rangedPotentialIterator.reset();
        maxPotentialRange = 0;
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (potentials[i].getRange() > maxPotentialRange) {
                    maxPotentialRange = potentials[i].getRange();
                }
            }
        }
    }

    /**
     * Adds the potential as a ranged potential that applies to the given 
     * AtomTypes.  This method creates a criterion for the potential and 
     * notifies the NeighborListManager of its existence.
     */
    protected void addRangedPotentialForTypes(IPotential potential, AtomType[] atomType) {
        super.addRangedPotentialForTypes(potential, atomType);
        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
    }

    public void removePotential(IPotential potential) {
        super.removePotential(potential);
        
        maxPotentialRange = 0;
        for (int i=0; i<allPotentials.length; i++) {
            double pRange = allPotentials[i].getRange();
            if (pRange == Double.POSITIVE_INFINITY) {
                continue;
            }
            if (pRange > maxPotentialRange) {
                maxPotentialRange = pRange;
            }
        }
    }

    private static final long serialVersionUID = 1L;
    private double range;
    private double maxPotentialRange;
}
