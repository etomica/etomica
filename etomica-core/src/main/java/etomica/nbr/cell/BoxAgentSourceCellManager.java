/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.box.Box;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * BoxAgentSource responsible for creating a NeighborCellManager.
 */
public class BoxAgentSourceCellManager implements BoxAgentSource<NeighborCellManager> {

    public BoxAgentSourceCellManager(IMoleculePositionDefinition positionDefinition, double range) {
        this.positionDefinition = positionDefinition;
        this.range = range;
    }
    
    public void setRange(double d) {
        range = d;
    }

    public NeighborCellManager makeAgent(Box box) {
        return new NeighborCellManager(box,range,positionDefinition);
    }
    
    public void releaseAgent(NeighborCellManager agent) {
    }
    
    protected double range;
    protected final IMoleculePositionDefinition positionDefinition;
}
