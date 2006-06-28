package etomica.nbr.cell;

import etomica.atom.Atom;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.SpeciesRoot;
import etomica.nbr.site.PotentialMasterSite;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.space.Space;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/*
 * Created on May 27, 2005
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

    /**
     * @param space
     * @param positionDefinition
     */
    public PotentialMasterCell(Space space, double range,
            AtomPositionDefinition positionDefinition) {
        this(space, range, new PhaseAgentSourceCellManager(positionDefinition));
    }
    
    public PotentialMasterCell(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource) {
        this(space, range, phaseAgentSource, new PhaseAgentManager(phaseAgentSource,null));
    }
    
    public PotentialMasterCell(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource,
            PhaseAgentManager agentManager) {
        super(space, phaseAgentSource, agentManager, new Api1ACell(space.D(),range,agentManager));
    }
    
    public double getRange() {
        return range;
    }
    
    public void setRange(double d) {
        ((Api1ACell)neighborIterator).getNbrCellIterator().setNeighborDistance(d);
        ((PhaseAgentSourceCellManager)phaseAgentSource).setRange(d);
        range = d;
    }
    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        PhaseAgentManager phaseAgentManager = getCellAgentManager();
        NeighborCellManager[] cellManagers = (NeighborCellManager[])phaseAgentManager.getAgents();
        NeighborCellManager manager = cellManagers[phase.getIndex()];
        manager.setPotentialRange(range);
        manager.setCellRange(getCellRange());
        return manager;
    }
    
    private double range;
}
