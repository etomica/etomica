package etomica.phase;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;

public class PhaseAgentSourceAtomManager implements PhaseAgentSource, java.io.Serializable {

    public PhaseAgentSourceAtomManager(AgentSource atomAgentSource) {
        super();
        this.atomAgentSource = atomAgentSource;
    }

    public Class getAgentClass() {
        return AtomAgentManager.class;
    }

    public Object makeAgent(Phase phase) {
        return new AtomAgentManager(atomAgentSource,phase);
    }

    public void releaseAgent(Object agent) {
        ((AtomAgentManager)agent).setPhase(null);
    }

    private final AgentSource atomAgentSource;
}
