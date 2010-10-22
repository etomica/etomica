package etomica.nbr.cell.molecule;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.IAtomPositionDefinition;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.ISpace;

/**
 * BoxAgentSource responsible for creating a NeighborCellManagerMolecular.
 *
 * @author Tai Boon Tan
 *
 */
public class BoxAgentSourceCellManagerMolecular implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceCellManagerMolecular(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        this.sim = sim;
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public void setRange(double d) {
        range = d;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerMolecular.class;
    }
    
    public Object makeAgent(IBox box) {
        NeighborCellManagerMolecular cellManager = new NeighborCellManagerMolecular(sim, box,range,positionDefinition, space);
        box.getBoundary().getEventManager().addListener(cellManager);
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