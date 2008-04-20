package etomica.box;

import etomica.api.IBox;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.BoxAgentManager.BoxAgentSource;

public class BoxAgentSourceAtomManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceAtomManager(AgentSource atomAgentSource) {
        super();
        this.atomAgentSource = atomAgentSource;
    }

    public Class getAgentClass() {
        return AtomLeafAgentManager.class;
    }

    public Object makeAgent(IBox box) {
        return new AtomLeafAgentManager(atomAgentSource,box);
    }

    public void releaseAgent(Object agent) {
        ((AtomLeafAgentManager)agent).dispose();
    }

    private static final long serialVersionUID = 1L;
    private final AgentSource atomAgentSource;
}
