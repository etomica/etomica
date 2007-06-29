/**
 * 
 */
package etomica.nbr.cell;

import etomica.atom.AtomPositionDefinition;
import etomica.box.Box;
import etomica.box.BoxAgentManager.BoxAgentSource;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManager(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    public void setRange(double d) {
        range = d;
    }
    
    public Class getAgentClass() {
        return NeighborCellManager.class;
    }
    
    public Object makeAgent(Box box) {
        NeighborCellManager cellManager = new NeighborCellManager(box,range,positionDefinition);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    private double range;
    private final AtomPositionDefinition positionDefinition;
}