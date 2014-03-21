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
