/**
 * 
 */
package etomica.threaded.domain;

import etomica.api.IBox;
import etomica.atom.AtomPositionDefinition;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.Space;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerThreaded implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManagerThreaded(AtomPositionDefinition positionDefinition, Space _space) {
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerThreaded.class;
    }
    
    public Object makeAgent(IBox box) {
        NeighborCellManagerThreaded cellManager = new NeighborCellManagerThreaded(box, 0, positionDefinition, space);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomPositionDefinition positionDefinition;
    private final Space space;
}