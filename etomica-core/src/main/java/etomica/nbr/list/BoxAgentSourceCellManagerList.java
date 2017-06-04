/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.IMoleculePositionDefinition;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.nbr.cell.BoxAgentSourceCellManager;
import etomica.nbr.cell.NeighborCellManager;
import etomica.space.Space;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerList extends BoxAgentSourceCellManager {

    public BoxAgentSourceCellManagerList(Simulation sim, IMoleculePositionDefinition positionDefinition, Space _space) {
        super(sim, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterList newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public NeighborCellManager makeAgent(Box box) {
        NeighborCellManagerList cellManager = new NeighborCellManagerList(sim, box,range,positionDefinition, space);
        cellManager.setPotentialMaster(potentialMaster);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }

    protected PotentialMasterList potentialMaster;
}
