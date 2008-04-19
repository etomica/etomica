package etomica.nbr.cell;

import etomica.api.IAtomPositionDefinition;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.ISpace;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManager(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        this.sim = sim;
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
        NeighborCellManager cellManager = new NeighborCellManager(sim, box,range,positionDefinition, space);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    protected final ISimulation sim;
    protected double range;
    protected final IAtomPositionDefinition positionDefinition;
    protected final ISpace space;
}