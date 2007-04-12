/**
 * 
 */
package etomica.threaded.domain;

import etomica.atom.AtomPositionDefinition;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;

/**
 * PhaseAgentSource responsible for creating a NeighborCellManager.
 */
public class PhaseAgentSourceCellManagerThreaded implements PhaseAgentSource, java.io.Serializable {

    public PhaseAgentSourceCellManagerThreaded(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerThreaded.class;
    }
    
    public Object makeAgent(Phase phase) {
        NeighborCellManagerThreaded cellManager = new NeighborCellManagerThreaded(phase, 0, positionDefinition);
        phase.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomPositionDefinition positionDefinition;
}