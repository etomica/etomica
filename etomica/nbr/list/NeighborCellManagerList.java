package etomica.nbr.list;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositionDefinition;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.nbr.cell.NeighborCellManager;
import etomica.space.ISpace;

/**
 * Subclass of NeighborCellManager that notifies the NeighborListManager when
 * an atom added to the box gets a cell.
 * 
 * @author Andrew Schultz
 */
public class NeighborCellManagerList extends NeighborCellManager {

    public NeighborCellManagerList(ISimulation sim, IBox box,
            double potentialRange, ISpace _space) {
        this(sim, box, potentialRange, null, _space);
    }

    public NeighborCellManagerList(ISimulation sim, IBox box,
            double potentialRange, IAtomPositionDefinition positionDefinition,
            ISpace _space) {
        super(sim, box, potentialRange, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterList newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public Object makeAgent(IAtomLeaf atom) {
        Object cell = super.makeAgent(atom);
        // oh, the humanity
        // we've made the cell and added the atom to it.  but the agent manager
        // won't actually have the cell until we return.  so explicitly set the
        // agent now.
        agentManager.setAgent(atom, cell);
        // now notify the NeighborListManager that the atom was added and has a
        // cell.  NeighborListManager wants to construct neighbor lists now.
        potentialMaster.getNeighborManager(box).addAtomNotify(atom);
        return cell;
    }

    protected PotentialMasterList potentialMaster;
}
