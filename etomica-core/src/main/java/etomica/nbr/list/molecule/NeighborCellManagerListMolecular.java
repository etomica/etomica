/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.atom.IMoleculePositionDefinition;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.simulation.Simulation;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.space.Space;

/**
 * Subclass of NeighborCellManager that notifies the NeighborListManager when
 * an atom added to the box gets a cell.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class NeighborCellManagerListMolecular extends NeighborCellManagerMolecular {

	public NeighborCellManagerListMolecular(Simulation sim, Box box,
                                            double potentialRange, Space _space) {
        this(sim, box, potentialRange, null, _space);
    }

    public NeighborCellManagerListMolecular(Simulation sim, Box box,
                                            double potentialRange, IMoleculePositionDefinition positionDefinition,
                                            Space _space) {
        super(sim, box, potentialRange, positionDefinition, _space);
    }

    public void setPotentialMaster(PotentialMasterListMolecular newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }

    public Object makeAgent(IMolecule molecule) {
        if (range == 0) {
            // no range, no lattice, etc
            return null;
        }
        Object cell = super.makeAgent(molecule);
        // if agentManager is null, we're in the constructor, which means we're
        // just starting up.  And we don't need to do this craziness.
        if (agentManager != null) {
            // oh, the humanity
            // we've made the cell and added the atom to it.  but the agent manager
            // won't actually have the cell until we return.  so explicitly set the
            // agent now.
            agentManager.setAgent(molecule, cell);
            // now notify the NeighborListManager that the atom was added and has a
            // cell.  NeighborListManager wants to construct neighbor lists now.
            potentialMaster.getNeighborManager(box).addMoleculeNotify(molecule);
        }
        return cell;
    }
    
    private static final long serialVersionUID = 1L;
    protected PotentialMasterListMolecular potentialMaster;
}
