package etomica.phase;

import java.lang.reflect.Array;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.simulation.Simulation;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;
import etomica.util.Arrays;

/**
 * AtomAgentManager acts on behalf of client classes (an AgentSource) to manage 
 * agents in every Atom in a phase.  When atoms are added or removed from the 
 * phase, the agents array (indexed by the atom's global index) is updated.  
 * The client should call getAgents() at any point where an atom might have 
 * have been added to the system because the old array would be stale at that
 * point. 
 * @author andrew
 */
public class PhaseAgentManager implements SimulationListener, java.io.Serializable {

    public PhaseAgentManager(PhaseAgentSource source, SpeciesRoot speciesRoot) {
        agentSource = source;
        setRoot(speciesRoot);
    }
    
    public Object[] getAgents() {
        return agents;
    }
    
    public void setRoot(SpeciesRoot newSpeciesRoot) {
        if (newSpeciesRoot == speciesRoot) {
            return;
        }
        if (speciesRoot != null) {
            // remove ourselves as a listener to the old phase
            speciesRoot.removeListener(this);
            for (int i=0; i<agents.length; i++) {
                if (agents[i] != null) {
                    agentSource.releaseAgent(agents[i]);
                }
            }
        }
        speciesRoot = newSpeciesRoot;
        if (speciesRoot == null) {
            // zero-out the array so the agents get GC'd
            agents = (Object[])Array.newInstance(agentSource.getAgentClass(),0);
            return;
        }
        speciesRoot.addListener(this);
        // hope the class returns an actual class with a null Atom and use it to construct
        // the array
        AtomIteratorArrayListSimple listIterator = new AtomIteratorArrayListSimple(((AtomTreeNodeGroup)speciesRoot.node).childList);
        listIterator.reset();
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),listIterator.size());
        while(listIterator.hasNext()) {
            addAgent(listIterator.nextAtom().node.parentPhase());
        }
    }
    
    public void actionPerformed(SimulationEvent evt) {
        if (evt.type() == SimulationEvent.PHASE_ADDED) {
            addAgent(evt.getPhase());
        }
        else if (evt.type() == SimulationEvent.PHASE_REMOVED) {
            Phase phase = evt.getPhase();
            int index = phase.getIndex()-1;
            agentSource.releaseAgent(agents[index]);
            for (int i=index; i<agents.length-1; i++) {
                agents[index] = agents[index+1];
            }
        }
    }
    
    protected void addAgent(Phase phase) {
        agents = Arrays.resizeArray(agents,phase.getIndex()+1);
        agents[phase.getIndex()] = agentSource.makeAgent(phase);
    }
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface PhaseAgentSource {
        public Class getAgentClass();
        
        public Object makeAgent(Phase phase);
        
        //allow any agent to be disconnected from other elements 
        public void releaseAgent(Object agent); 
    }

    private final PhaseAgentSource agentSource;
    protected Object[] agents;
    private SpeciesRoot speciesRoot;
}
