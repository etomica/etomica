/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.simulation.Simulation;

/**
 * @author Tai Boon Tan
 */
public class BoxAgentSourceMoleculeManager implements BoxAgentSource<MoleculeAgentManager> {

    private final MoleculeAgentSource moleculeAgentSource;
    protected Simulation sim;

    public BoxAgentSourceMoleculeManager(MoleculeAgentSource moleculeAgentSource, Simulation sim) {
        super();
        this.sim = sim;
        this.moleculeAgentSource = moleculeAgentSource;
    }

    public MoleculeAgentManager makeAgent(Box box) {
        return new MoleculeAgentManager(sim, box, moleculeAgentSource);
    }

    public void releaseAgent(MoleculeAgentManager agent) {
        agent.dispose();
    }
}
