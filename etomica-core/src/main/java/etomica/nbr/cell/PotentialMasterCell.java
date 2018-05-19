/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.site.PotentialMasterSite;
import etomica.simulation.Simulation;
import etomica.space.Space;

import java.util.function.Function;

/**
 * A PotentialMaster for use with a Simulation where cell-listing of atoms and
 * their neighbors is appropriate.
 *
 * @author Andrew Schultz
 */
public class PotentialMasterCell extends PotentialMasterSite {

    private double range;
    private final BoxAgentSourceCellManager boxAgentSource;
    private final BoxAgentManager<NeighborCellManager> boxAgentManagerNeighborCell;

    /**
     * Creates PotentialMasterCell with default (1.0) range.  Range
     * should be set manually via setRange method.
     */
    public PotentialMasterCell(Simulation sim, Space _space) {
        this(sim, 1.0, _space);
    }

    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     *
     * @param _space the governing Space
     * @param range  the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(Simulation sim, double range, Space _space) {
        this(sim, range, (IMoleculePositionDefinition) null, _space);
    }

    public PotentialMasterCell(Simulation sim, double range,
                               IMoleculePositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManager(sim, positionDefinition, _space), _space);
    }

    public PotentialMasterCell(Simulation sim, double range,
                               BoxAgentSourceCellManager boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<>(boxAgentSource, sim), _space);
    }

    public PotentialMasterCell(Simulation sim, double range, BoxAgentSourceCellManager boxAgentSource,
                               BoxAgentManager<NeighborCellManager> agentManager, Space _space) {
        super(sim, agentManager, new Api1ACell(_space.D(), range, agentManager));
        this.boxAgentSource = boxAgentSource;
        this.boxAgentManagerNeighborCell = agentManager;
        setRange(range);
    }

    public double getRange() {
        return range;
    }

    public void setRange(double d) {
        ((Api1ACell) neighborIterator).getNbrCellIterator().setNeighborDistance(d);
        boxAgentSource.setRange(d);
        range = d;

        for (NeighborCellManager neighborCellManager : boxAgentManagerNeighborCell.getAgents().values()) {
            neighborCellManager.setPotentialRange(range);
        }
    }

    public void setCellRange(int d) {
        super.setCellRange(d);

        for (NeighborCellManager neighborCellManager : this.boxAgentManagerNeighborCell.getAgents().values()) {
            neighborCellManager.setCellRange(this.getCellRange());
        }
    }

    public NeighborCellManager getNbrCellManager(Box box) {
        NeighborCellManager manager = boxAgentManagerNeighborCell.getAgent(box);
        manager.setPotentialRange(range);
        int cr = getCellRange();
        if (cr == 0) throw new RuntimeException("need to set cell range first");
        manager.setCellRange(cr);
        return manager;
    }

    /**
     * Reassign atoms to cell lists for all boxes.
     */
    public void reset() {
        boxAgentManager.getAgents().values().forEach(BoxCellManager::assignCellAll);
    }
}
