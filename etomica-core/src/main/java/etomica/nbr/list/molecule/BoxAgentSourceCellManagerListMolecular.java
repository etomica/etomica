/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.box.Box;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.cell.molecule.BoxAgentSourceCellManagerMolecular;
import etomica.space.Space;
import etomica.species.SpeciesManager;

/**
 * BoxAgentSource responsible for creating a NeighborCellManagerMolecular.
 *
 * @author Tai Boon Tan
 */
public class BoxAgentSourceCellManagerListMolecular extends BoxAgentSourceCellManagerMolecular {

    public BoxAgentSourceCellManagerListMolecular(SpeciesManager sm, IMoleculePositionDefinition positionDefinition, Space _space) {
        super(sm, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterListMolecular newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public NeighborCellManagerListMolecular makeAgent(Box box) {
        NeighborCellManagerListMolecular cellManager = new NeighborCellManagerListMolecular(sm, box, range, positionDefinition, space);
        cellManager.setPotentialMaster(potentialMaster);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }

    protected PotentialMasterListMolecular potentialMaster;
}
