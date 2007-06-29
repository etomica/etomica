/**
 * 
 */
package etomica.threaded.domain;

import etomica.atom.AtomPositionDefinition;
import etomica.box.Box;
import etomica.box.BoxAgentManager.BoxAgentSource;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerThreaded implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManagerThreaded(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerThreaded.class;
    }
    
    public Object makeAgent(Box box) {
        NeighborCellManagerThreaded cellManager = new NeighborCellManagerThreaded(box, 0, positionDefinition);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomPositionDefinition positionDefinition;
}