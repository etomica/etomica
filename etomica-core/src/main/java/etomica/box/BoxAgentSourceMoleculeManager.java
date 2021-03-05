/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.simulation.Simulation;

/**
 * BoxAgentSource that creates an MoleculeAgentManager.
 * @author Tai Boon Tan
 */
public class BoxAgentSourceMoleculeManager<E> implements BoxAgentSource<MoleculeAgentManager<E>> {

    private final MoleculeAgentSource<E> moleculeAgentSource;
    protected Simulation sim;

    /**
     * @param moleculeAgentSource object that makes the molecule agents
     * @param sim the governing simulation
     */
    public BoxAgentSourceMoleculeManager(MoleculeAgentSource<E> moleculeAgentSource, Simulation sim) {
        super();
        this.sim = sim;
        this.moleculeAgentSource = moleculeAgentSource;
    }

    public MoleculeAgentManager<E> makeAgent(Box box) {
        return new MoleculeAgentManager<>(sim.getSpeciesManager(), box, moleculeAgentSource);
    }

    public void releaseAgent(MoleculeAgentManager<E> agent) {
        agent.dispose();
    }
}
