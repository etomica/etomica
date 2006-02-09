/**
 * 
 */
package etomica.nbr.cell;

import etomica.atom.AtomPositionDefinition;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;

public class PhaseAgentSourceCellManager implements PhaseAgentSource, java.io.Serializable {
    public PhaseAgentSourceCellManager(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    public void setRange(double d) {
        range = d;
    }
    
    public Class getAgentClass() {
        return NeighborCellManager.class;
    }
    
    public Object makeAgent(Phase phase) {
        return new NeighborCellManager(phase,range,positionDefinition);
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private double range;
    private final AtomPositionDefinition positionDefinition;
}