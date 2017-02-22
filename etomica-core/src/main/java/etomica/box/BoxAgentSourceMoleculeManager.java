/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.box.BoxAgentManager.BoxAgentSource;

/**
 * @author Tai Boon Tan
 *
 */
public class BoxAgentSourceMoleculeManager implements BoxAgentSource<MoleculeAgentManager> {

    public BoxAgentSourceMoleculeManager(MoleculeAgentSource moleculeAgentSource, ISimulation sim) {
        super();
        this.sim = sim;
        this.moleculeAgentSource = moleculeAgentSource;
    }

    public MoleculeAgentManager makeAgent(IBox box) {
        return new MoleculeAgentManager(sim, box, moleculeAgentSource);
    }

    public void releaseAgent(MoleculeAgentManager agent) {
        agent.dispose();
    }

    protected ISimulation sim;
    private final MoleculeAgentSource moleculeAgentSource;
}
