package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.util.IRandom;

public class ColorSchemeRandom extends ColorScheme implements AgentSource {
    
    public ColorSchemeRandom(Phase phase, IRandom random) {
        this.random = random;
        agentManager = new AtomAgentManager(this, phase);
    }
    
    public Color getAtomColor(AtomLeaf a) {
        return (Color)agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(Atom a) {
        return new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
    }
    
    public void releaseAgent(Object agent, Atom atom) {}

    private static final long serialVersionUID = 2L;
    private final AtomAgentManager agentManager;
    private final IRandom random;
}