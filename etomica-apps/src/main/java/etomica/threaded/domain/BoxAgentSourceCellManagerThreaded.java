/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.domain;

import etomica.atom.IMoleculePositionDefinition;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.space.Space;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManagerThreaded implements BoxAgentSource<NeighborCellManagerThreaded> {

    public BoxAgentSourceCellManagerThreaded(Simulation sim, IMoleculePositionDefinition positionDefinition, Space _space) {
        this.sim = sim;
        this.positionDefinition = positionDefinition;
        this.space = _space;
    }
    
    public Class getAgentClass() {
        return NeighborCellManagerThreaded.class;
    }
    
    public NeighborCellManagerThreaded makeAgent(Box box) {
        NeighborCellManagerThreaded cellManager = new NeighborCellManagerThreaded(sim, box, 0, positionDefinition, space);
        box.getBoundary().getEventManager().addListener(cellManager);
        return cellManager;
    }
    
    public void releaseAgent(NeighborCellManagerThreaded agent) {
    }
    
    protected final Simulation sim;
    private final IMoleculePositionDefinition positionDefinition;
    private final Space space;
}
