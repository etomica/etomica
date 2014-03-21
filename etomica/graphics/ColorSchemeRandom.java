package etomica.graphics;
import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;

public class ColorSchemeRandom extends ColorScheme implements AgentSource<Color> {
    
    public ColorSchemeRandom(IBox box, IRandom random) {
    	super();
        this.random = random;
        agentManager = new AtomLeafAgentManager<Color>(this, box, Color.class);
    }
    
    public Color getAtomColor(IAtom a) {
        return agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Color makeAgent(IAtom a) {
        return new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
    }
    
    public void releaseAgent(Color agent, IAtom atom) {}

    private static final long serialVersionUID = 2L;
    private final AtomLeafAgentManager<Color> agentManager;
    private final IRandom random;
}