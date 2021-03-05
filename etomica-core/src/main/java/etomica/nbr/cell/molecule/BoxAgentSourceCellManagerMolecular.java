/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell.molecule;

import etomica.box.Box;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.space.Space;
import etomica.species.SpeciesManager;

/**
 * BoxAgentSource responsible for creating a NeighborCellManagerMolecular.
 *
 * @author Tai Boon Tan
 */
public class BoxAgentSourceCellManagerMolecular implements BoxAgentSource<NeighborCellManagerMolecular> {

    public BoxAgentSourceCellManagerMolecular(SpeciesManager sm, IMoleculePositionDefinition positionDefinition, Space _space) {
        this.sm = sm;
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }

    public void setRange(double d) {
        range = d;
    }

    public NeighborCellManagerMolecular makeAgent(Box box) {
        NeighborCellManagerMolecular cellManager = new NeighborCellManagerMolecular(sm, box, range, positionDefinition, space);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }

    public void releaseAgent(NeighborCellManagerMolecular agent) {
    }

    protected final SpeciesManager sm;
    protected double range;
    protected final IMoleculePositionDefinition positionDefinition;
    protected final Space space;
}
