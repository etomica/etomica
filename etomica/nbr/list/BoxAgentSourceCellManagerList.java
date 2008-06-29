package etomica.nbr.list;

import etomica.api.IAtomPositionDefinition;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.nbr.cell.BoxAgentSourceCellManager;
import etomica.space.ISpace;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerList extends BoxAgentSourceCellManager {

    public BoxAgentSourceCellManagerList(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        super(sim, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterList newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public Object makeAgent(IBox box) {
        NeighborCellManagerList cellManager = new NeighborCellManagerList(sim, box,range,positionDefinition, space);
        cellManager.setPotentialMaster(potentialMaster);
        box.getEventManager().addListener(cellManager);
        return cellManager;
    }

    private static final long serialVersionUID = 1L;
    protected PotentialMasterList potentialMaster;
}