package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;

/**
 * Parent class for color schemes that are best implemented by attaching colors
 * to all atoms at once, rather than determining the color of each as it is drawn.
 * The colorAllAtoms method is called by the display if it determines that the
 * ColorScheme is a subclass of this one.
 */
public abstract class ColorSchemeCollective extends ColorScheme implements AgentSource {
    
    protected AtomAgentManagerColor agentManager;
    
    public ColorSchemeCollective(Phase phase) {
        agentManager = new AtomAgentManagerColor(this, phase);
    }
    
    //determine color
    //then assign it to atom like this: atomColors[atomIndex] = color
    public abstract void colorAllAtoms();
    
    public Color getAtomColor(AtomLeaf a) {
        return (Color)agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(Atom a) {
        return null;
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
    
    /**
     * Inner class to let us cheat and access and modify elements of the agents array.
     */
    protected static class AtomAgentManagerColor extends AtomAgentManager {
        public AtomAgentManagerColor(ColorSchemeCollective colorScheme, Phase phase) {
            super(colorScheme, phase, false);
        }

        /**
         * Returns the array of Colors for Atoms in the Phase, indexed by the Atom's
         * global index.
         */
        protected Color[] getAgents() {
            return (Color[])agents;
        }
    }
}
