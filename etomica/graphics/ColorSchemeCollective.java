package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.SpeciesRoot;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;
import etomica.util.Arrays;

/**
 * Parent class for color schemes that are best implemented by attaching colors
 * to all atoms at once, rather than determining the color of each as it is drawn.
 * The colorAllAtoms method is called by the display if it determines that the
 * ColorScheme is a subclass of this one.
 */
    
public abstract class ColorSchemeCollective extends ColorScheme implements AgentSource, SimulationListener {
    
    private AtomAgentManager[] agentManager;
    protected Color[] atomColors;
    
    public ColorSchemeCollective() {
        super();
    }
    
    //determine color
    //then assign it to atom like this: atom.allatomAgents[agentIndex] = color
    public void colorAllAtoms(Phase phase){
        if (agentManager == null) {
            ((SpeciesRoot)phase.getSpeciesMaster().node.parentGroup()).addListener(this);
            agentManager = new AtomAgentManager[0];
        }
        int index = phase.getIndex();
        if (index+1 > agentManager.length) {
            agentManager = (AtomAgentManager[])Arrays.resizeArray(agentManager,index+1);
        }
        if (agentManager[index] == null) {
            agentManager[index] = new AtomAgentManager(this,phase);
        }
        atomColors = (Color[])agentManager[index].getAgents();
    }
    
    public Color getAtomColor(AtomLeaf a) {return atomColors[a.getGlobalIndex()];}
   
    public void actionPerformed(SimulationEvent evt) {
        if (evt.type() == SimulationEvent.PHASE_REMOVED) {
            Phase p = evt.getPhase();
            if (agentManager[p.getIndex()] != null) {
                // allow AtomAgentManager to register itself as a PhaseListener
                agentManager[p.getIndex()].setPhase(null);
                agentManager[p.getIndex()] = null;
            }
            // compact the array if there are null elements at the end.
            int i;
            for (i=agentManager.length; i>0; i++) {
                if (agentManager[i-1] != null) {
                    break;
                }
            }
            if (i > 0 && i < agentManager.length) {
                agentManager = (AtomAgentManager[])Arrays.resizeArray(agentManager,i);
            }
        }
    }
    
    //set aside a agent index entry to store the color with the atom 
    public Object makeAgent(Atom a) {
        // return a color if a is null (AtomAgentManager uses this to determine
        // the class
        if (a == null) {return Color.white;}
        return null;
    }
    
    public void releaseAgent(Object agent) {}
}

