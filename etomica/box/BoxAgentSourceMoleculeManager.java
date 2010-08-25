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
public class BoxAgentSourceMoleculeManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceMoleculeManager(MoleculeAgentSource moleculeAgentSource, ISimulation sim) {
        super();
        this.sim = sim;
        this.moleculeAgentSource = moleculeAgentSource;
    }

    public Class getAgentClass() {
        return MoleculeAgentManager.class;
    }

    public Object makeAgent(IBox box) {
        return new MoleculeAgentManager(sim, box, moleculeAgentSource);
    }

    public void releaseAgent(Object agent) {
        ((MoleculeAgentManager)agent).dispose();
    }

    protected ISimulation sim;
    private static final long serialVersionUID = 1L;
    private final MoleculeAgentSource moleculeAgentSource;
}
