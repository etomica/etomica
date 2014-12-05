/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.IAtomPositionDefinition;
import etomica.nbr.cell.molecule.BoxAgentSourceCellManagerMolecular;
import etomica.space.ISpace;

/**
 * BoxAgentSource responsible for creating a NeighborCellManagerMolecular.
 * 
 * @author Tai Boon Tan
 *
 */
public class BoxAgentSourceCellManagerListMolecular extends BoxAgentSourceCellManagerMolecular {

    public BoxAgentSourceCellManagerListMolecular(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        super(sim, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterListMolecular newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public NeighborCellManagerListMolecular makeAgent(IBox box) {
        NeighborCellManagerListMolecular cellManager = new NeighborCellManagerListMolecular(sim, box,range,positionDefinition, space);
        cellManager.setPotentialMaster(potentialMaster);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }

    protected PotentialMasterListMolecular potentialMaster;
}