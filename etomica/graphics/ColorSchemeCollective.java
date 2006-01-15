package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.simulation.Simulation;

/**
 * Parent class for color schemes that are best implemented by attaching colors
 * to all atoms at once, rather than determining the color of each as it is drawn.
 * The colorAllAtoms method is called by the display if it determines that the
 * ColorScheme is a subclass of this one.
 */
    
public abstract class ColorSchemeCollective extends ColorScheme implements AgentSource {
    
    private AtomAgentManager[] agentManager;
    protected Color[] atomColors;
    private final PhaseAgentManager phaseAgentManager;
    
    public ColorSchemeCollective(Simulation sim) {
        super();
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.speciesRoot);
    }
    
    //determine color
    //then assign it to atom like this: atomColors[atomIndex] = color
    public void colorAllAtoms(Phase phase){
        agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        atomColors = (Color[])agentManager[phase.getIndex()].getAgents();
    }
    
    public Color getAtomColor(AtomLeaf a) {return atomColors[a.getGlobalIndex()];}
   
    public Object makeAgent(Atom a) {
        // return a color if a is null (AtomAgentManager uses this to determine
        // the class)
        if (a == null) {return Color.white;}
        return null;
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
}

