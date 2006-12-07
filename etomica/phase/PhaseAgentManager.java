package etomica.phase;

import java.lang.reflect.Array;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;
import etomica.simulation.SimulationPhaseAddedEvent;
import etomica.simulation.SimulationPhaseEvent;
import etomica.util.Arrays;

/**
 * PhaseAgentManager acts on behalf of client classes (a PhaseAgentSource) to manage 
 * agents for each Phase in a simulation.  When Phase instances are added or removed 
 * from the simulation, the agents array (indexed by the phase's index) is updated.  
 * The client should call getAgents() at any point where a Phase might have have been 
 * added to (or removed from) the system because the old array would be stale at that
 * point. 
 * @author andrew
 */
public class PhaseAgentManager implements SimulationListener, java.io.Serializable {

    public PhaseAgentManager(PhaseAgentSource source, SpeciesRoot speciesRoot) {
        this(source, speciesRoot, true);
    }

    public PhaseAgentManager(PhaseAgentSource source, SpeciesRoot speciesRoot, boolean isBackend) {
        agentSource = source;
        if (speciesRoot != null) {
            setRoot(speciesRoot);
        }
        this.isBackend = isBackend;
    }
    
    /**
     * Returns the agent associated with the given phase
     */
    public Object getAgent(Phase phase) {
        return agents[phase.getIndex()];
    }
    
    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Sets the species root parent in the Simulation for Phases to be tracked.
     */
    public void setRoot(SpeciesRoot newSpeciesRoot) {
        if (newSpeciesRoot == speciesRoot) {
            return;
        }
        speciesRoot = newSpeciesRoot;
        speciesRoot.getEventManager().addListener(this, isBackend);
        // hope the class returns an actual class with a null Atom and use it to construct
        // the array
        AtomIteratorArrayListSimple listIterator = new AtomIteratorArrayListSimple(((AtomTreeNodeGroup)speciesRoot.getNode()).childList);
        listIterator.reset();
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),listIterator.size());
        while(listIterator.hasNext()) {
            addAgent(listIterator.nextAtom().getNode().parentPhase());
        }
    }
    
    /**
     * Notifies the PhaseAgentManager that it should release all agents and 
     * stop listening for events from the simulation.
     */
    public void dispose() {
        // remove ourselves as a listener to the old phase
        speciesRoot.getEventManager().removeListener(this);
        for (int i=0; i<agents.length; i++) {
            if (agents[i] != null) {
                agentSource.releaseAgent(agents[i]);
            }
        }
        agents = null;
    }
    
    public void actionPerformed(SimulationEvent evt) {
        if (evt instanceof SimulationPhaseAddedEvent) {
            addAgent(((SimulationPhaseEvent)evt).getPhase());
        }
        else if (evt instanceof SimulationPhaseAddedEvent) {
            Phase phase = ((SimulationPhaseEvent)evt).getPhase();
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

    private static final long serialVersionUID = 1L;
    private final PhaseAgentSource agentSource;
    protected Object[] agents;
    private SpeciesRoot speciesRoot;
    private final boolean isBackend;
    
    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator {
        protected AgentIterator(PhaseAgentManager agentManager) {
            this.agentManager = agentManager;
        }
        
        public void reset() {
            cursor = 0;
            agents = agentManager.agents;
        }
        
        public boolean hasNext() {
            while (cursor < agents.length) {
                if (agents[cursor] != null) {
                    return true;
                }
                cursor++;
            }
            return false;
        }
        
        public Object next() {
            cursor++;
            while (cursor-1 < agents.length) {
                if (agents[cursor-1] != null) {
                    return agents[cursor-1];
                }
                cursor++;
            }
            return null;
        }
        
        private final PhaseAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
}
