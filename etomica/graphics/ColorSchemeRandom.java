package etomica.graphics;
import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;

public class ColorSchemeRandom extends ColorScheme implements AgentSource {
    
    public ColorSchemeRandom(IBox box, IRandom random) {
        this.random = random;
        agentManager = new AtomLeafAgentManager(this, box);
    }
    
    public Color getAtomColor(IAtom a) {
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