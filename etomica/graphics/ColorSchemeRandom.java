package etomica.graphics;
import java.awt.Color;

import etomica.api.IAtomLeaf;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;

public class ColorSchemeRandom extends ColorScheme implements AgentSource {
    
    public ColorSchemeRandom(ISimulation sim, IBox box, IRandom random) {
    	super(sim);
        this.random = random;
        agentManager = new AtomLeafAgentManager(this, box);
    }
    
    public Color getAtomColor(IAtomLeaf a) {
        return (Color)agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(IAtomLeaf a) {
        return new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}

    private static final long serialVersionUID = 2L;
    private final AtomLeafAgentManager agentManager;
    private final IRandom random;
}