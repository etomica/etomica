/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.atom.IMoleculePositionDefinition;
import etomica.nbr.cell.Cell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.space.Space;

/**
 * Subclass of NeighborCellManager that notifies the NeighborListManager when
 * an atom added to the box gets a cell.
 * 
 * @author Andrew Schultz
 */
public class NeighborCellManagerList extends NeighborCellManager {

    public NeighborCellManagerList(Simulation sim, Box box,
                                   double potentialRange, Space _space) {
        this(sim, box, potentialRange, null, _space);
    }

    public NeighborCellManagerList(Simulation sim, Box box,
                                   double potentialRange, IMoleculePositionDefinition positionDefinition,
                                   Space _space) {
        super(sim, box, potentialRange, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterList newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public Cell makeAgent(IAtom atom, Box agentBox) {
        if (range == 0) {
            // no range, no lattice, etc
            return null;
        }
        Cell cell = super.makeAgent(atom, agentBox);
        // if agentManager is null, we're in the constructor, which means we're
        // just starting up.  And we don't need to do this craziness.
        if (agentManager != null) {
            // oh, the humanity
            // we've made the cell and added the atom to it.  but the agent manager
            // won't actually have the cell until we return.  so explicitly set the
            // agent now.
            agentManager.setAgent(atom, cell);
            // now notify the NeighborListManager that the atom was added and has a
            // cell.  NeighborListManager wants to construct neighbor lists now.
            potentialMaster.getNeighborManager(box).addAtomNotify(atom);
        }
        return cell;
    }

    protected PotentialMasterList potentialMaster;
}
