/**
 * 
 */
package etomica.nbr.cell;

import etomica.api.IBox;
import etomica.atom.AtomPositionDefinition;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.Space;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManager(AtomPositionDefinition positionDefinition, Space _space) {
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public void setRange(double d) {
        range = d;
    }
    
    public Class getAgentClass() {
        return NeighborCellManager.class;
    }
    
    public Object makeAgent(IBox box) {
        NeighborCellManager cellManager = new NeighborCellManager(box,range,positionDefinition, space);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    private double range;
    private final AtomPositionDefinition positionDefinition;
    private final Space space;
}