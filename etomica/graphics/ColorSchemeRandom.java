package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;

public class ColorSchemeRandom extends ColorScheme implements AgentSource {
    
    public ColorSchemeRandom(Phase phase) {
        agentManager = new AtomAgentManager(this, phase);
    }
    
    public Color getAtomColor(AtomLeaf a) {
        return (Color)agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(Atom a) {
        return ConstantsGraphic.randomColor();
    }
    
    public void releaseAgent(Object agent, Atom atom) {}

    private static final long serialVersionUID = 1L;
    private final AtomAgentManager agentManager;
}