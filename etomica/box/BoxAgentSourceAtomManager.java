package etomica.box;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.BoxAgentManager.BoxAgentSource;

public class BoxAgentSourceAtomManager implements BoxAgentSource, java.io.Serializable {

    public BoxAgentSourceAtomManager(AgentSource atomAgentSource) {
        super();
        this.atomAgentSource = atomAgentSource;
    }

    public Class getAgentClass() {
        return AtomAgentManager.class;
    }

    public Object makeAgent(Box box) {
        return new AtomAgentManager(atomAgentSource,box);
    }

    public void releaseAgent(Object agent) {
        ((AtomAgentManager)agent).dispose();
    }

    private static final long serialVersionUID = 1L;
    private final AgentSource atomAgentSource;
}
