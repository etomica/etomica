/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.IAtomPositionDefinition;
import etomica.nbr.cell.BoxAgentSourceCellManager;
import etomica.nbr.cell.NeighborCellManager;
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

    public NeighborCellManager makeAgent(IBox box) {
        NeighborCellManagerList cellManager = new NeighborCellManagerList(sim, box,range,positionDefinition, space);
        cellManager.setPotentialMaster(potentialMaster);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }

    protected PotentialMasterList potentialMaster;
}