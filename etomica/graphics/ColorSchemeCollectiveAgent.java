package etomica.graphics;
import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;

/**
 * Parent class for color schemes that are best implemented by attaching colors
 * to all atoms at once, rather than determining the color of each as it is drawn.
 * The colorAllAtoms method is called by the display if it determines that the
 * ColorScheme is a subclass of this one.
 */
public abstract class ColorSchemeCollectiveAgent extends ColorScheme implements AgentSource, ColorSchemeCollective {
    
    protected AtomLeafAgentManager agentManager;
    
    public ColorSchemeCollectiveAgent(ISimulation sim, IBox box) {
    	super(sim);
        agentManager = new AtomLeafAgentManager(this, box);
    }
    
    //determine color
    //then assign it to atom like this: atomColors[atomIndex] = color
    /* (non-Javadoc)
     * @see etomica.graphics.IColorSchemeCollective#colorAllAtoms()
     */
    public abstract void colorAllAtoms();
    
    public Color getAtomColor(IAtom a) {
        return (Color)agentManager.getAgent((IAtomLeaf)a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(IAtomLeaf a) {
        return null;
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}
    
}
