package etomica.graphics;
import java.awt.Color;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.util.IRandom;

public class ColorSchemeRandom extends ColorScheme implements AgentSource {
    
    public ColorSchemeRandom(Phase phase, IRandom random) {
        this.random = random;
        agentManager = new AtomLeafAgentManager(this, phase);
    }
    
    public Color getAtomColor(AtomLeaf a) {
        return (Color)agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(IAtom a) {
        return new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}

    private static final long serialVersionUID = 2L;
    private final AtomLeafAgentManager agentManager;
    private final IRandom random;
}