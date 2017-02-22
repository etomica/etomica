/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
public class BoxAgentSourceCellManagerMolecular implements BoxAgentSource<NeighborCellManagerMolecular> {

    public BoxAgentSourceCellManagerMolecular(ISimulation sim, IAtomPositionDefinition positionDefinition, ISpace _space) {
        this.sim = sim;
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public void setRange(double d) {
        range = d;
    }

    public NeighborCellManagerMolecular makeAgent(IBox box) {
        NeighborCellManagerMolecular cellManager = new NeighborCellManagerMolecular(sim, box,range,positionDefinition, space);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(NeighborCellManagerMolecular agent) {
    }
    
    protected final ISimulation sim;
    protected double range;
    protected final IAtomPositionDefinition positionDefinition;
    protected final ISpace space;
}