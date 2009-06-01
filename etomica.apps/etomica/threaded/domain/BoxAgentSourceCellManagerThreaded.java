package etomica.threaded.domain;

import etomica.api.IAtomPositionDefinition;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.ISpace;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerThreaded implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManagerThreaded(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        this.sim = sim;
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerThreaded.class;
    }
    
    public Object makeAgent(IBox box) {
        NeighborCellManagerThreaded cellManager = new NeighborCellManagerThreaded(sim, box, 0, positionDefinition, space);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(Object agent) {
    }
    
    private static final long serialVersionUID = 1L;
    protected final ISimulation sim;
    private final IAtomPositionDefinition positionDefinition;
    private final ISpace space;
}