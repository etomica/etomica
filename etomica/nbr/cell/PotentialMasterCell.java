package etomica.nbr.cell;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomType;
import etomica.nbr.site.PotentialMasterSite;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
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
    public PotentialMasterCell(Space space) {
        this(space,1.0);
    }
    
    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     * 
     * @param space the governing Space
     * @param range the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(Space space, double range) {
        this(space, range, (AtomPositionDefinition)null);
    }

    public PotentialMasterCell(Space space, double range,
            AtomPositionDefinition positionDefinition) {
        this(space, range, new PhaseAgentSourceCellManager(positionDefinition));
    }
    
    public PotentialMasterCell(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource) {
        this(space, range, phaseAgentSource, new PhaseAgentManager(phaseAgentSource));
    }
    
    public PotentialMasterCell(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource,
            PhaseAgentManager agentManager) {
        super(space, phaseAgentSource, agentManager, new Api1ACell(space.D(),range,agentManager));
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
    
    public void setRange(double d) {
        ((Api1ACell)neighborIterator).getNbrCellIterator().setNeighborDistance(d);
        ((PhaseAgentSourceCellManager)phaseAgentSource).setRange(d);
        range = d;
    }
    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        NeighborCellManager manager = (NeighborCellManager)phaseAgentManager.getAgent(phase);
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
            Potential[] potentials = potentialArray.getPotentials();
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
    protected void addRangedPotentialForTypes(Potential potential, AtomType[] atomType) {
        super.addRangedPotentialForTypes(potential, atomType);
        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
    }

    public void removePotential(Potential potential) {
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
