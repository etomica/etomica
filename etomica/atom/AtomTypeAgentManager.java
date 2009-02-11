package etomica.atom;

import java.io.Serializable;
import java.lang.reflect.Array;

import etomica.api.IAtomTypeLeaf;
import etomica.api.IEvent;
import etomica.api.IEventManager;
import etomica.api.IListener;
import etomica.api.ISpecies;
import etomica.api.ISpeciesManager;
import etomica.simulation.SimulationAtomTypeIndexChangedEvent;
import etomica.simulation.SimulationAtomTypeMaxIndexEvent;
import etomica.simulation.SimulationSpeciesAddedEvent;
import etomica.simulation.SimulationSpeciesRemovedEvent;
import etomica.util.Arrays;

/**
 * AtomTypeAgentManager acts on behalf of client classes (an AgentSource) to 
 * manage agents in every AtomType in a box.  When species are added or 
 * removed from the simulation, the agents array (indexed by the AtomType's 
 * global index) is updated.  The client should call getAgents() at any point 
 * where an atom might have have been added to the system because the old array
 * would be stale at that point.
 * @author andrew
 */
public class AtomTypeAgentManager implements IListener, java.io.Serializable {

    public AtomTypeAgentManager(AgentSource source) {
        agentSource = source;
        isBackend = true;
    }
    
    public AtomTypeAgentManager(AgentSource source, ISpeciesManager speciesManager,
            IEventManager simEventManager, boolean isBackend) {
        agentSource = source;
        this.isBackend = isBackend;
        init(speciesManager, simEventManager);
    }        
    
    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Sets the agent associated with the given atom type to be the given
     * agent.  The AtomType must be from the ISimulation.  The AtomType's old
     * agent is not "released".  This should be done manually if needed.
     */
    public void setAgent(IAtomTypeLeaf atomType, Object newAgent) {
        agents[atomType.getIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the agent the given AtomType.  For repeated
     * access to the agents from multiple AtomTypes, it might be faster to use 
     * the above getAgents method.
     */
    public Object getAgent(IAtomTypeLeaf type) {
        return agents[type.getIndex()];
    }
    
    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(ISpecies parentType) {
        for (int i=0; i<parentType.getChildTypeCount(); i++) {
            IAtomTypeLeaf leafType = parentType.getChildType(i);
            Object agent = agents[leafType.getIndex()];
            if (agent != null) {
                agentSource.releaseAgent(agent, leafType);
                agents[leafType.getIndex()] = null;
            }
        }
    }
    
    private void makeAllAgents() {
        for (int i=0; i<speciesManager.getSpeciesCount(); i++) {
            ISpecies parentType = speciesManager.getSpecies(i);
            for (int j=0; j<parentType.getChildTypeCount(); j++) {
                addAgent(parentType.getChildType(j));
            }
        }
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private int getGlobalMaxIndex() {
        int max = 0;
        for (int i=0; i<speciesManager.getSpeciesCount(); i++) {
            if (speciesManager.getSpecies(i).getIndex() > max) {
                max = speciesManager.getSpecies(i).getIndex();
            }
            int childMax = getMaxIndexOfChildren(speciesManager.getSpecies(i));
            if (childMax > max) {
                max = childMax;
            }
        }
        return max;
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private static int getMaxIndexOfChildren(ISpecies parentType) {
        int max = 0;
        for (int i=0; i<parentType.getChildTypeCount(); i++) {
            if (parentType.getChildType(i).getIndex() > max) {
                max = parentType.getChildType(i).getIndex();
            }
        }
        return max;
    }
    
    /**
     * Unregisters this class as a listener for AtomType-related events and 
     * releases its agents.
     */
    public void dispose() {
        // remove ourselves as a listener to the old box
        simEventManager.removeListener(this);
        for (int i=0; i<speciesManager.getSpeciesCount(); i++) {
            releaseAgents(speciesManager.getSpecies(i));
        }
        agents = null;
    }
    
    /**
     * Sets the SpeciesRoot for which this AtomAgentManager will manage 
     * AtomType agents.
     */
    public void init(ISpeciesManager newSpeciesManager, IEventManager newSimEventManager) {
        simEventManager = newSimEventManager;
        speciesManager = newSpeciesManager;
        simEventManager.addListener(this, isBackend);

        int numTypes = getGlobalMaxIndex()+1;
        
        agents = (Object[])Array.newInstance(agentSource.getSpeciesAgentClass(), numTypes);
        // fill in the array with agents from all the atoms
        makeAllAgents();
    }
    
    public void actionPerformed(IEvent evt) {
        // we learn about new Species via AtomTypeAdded events
        if (evt instanceof SimulationSpeciesRemovedEvent) {
            releaseAgents(((SimulationSpeciesRemovedEvent)evt).getSpecies());
        }
        else if (evt instanceof SimulationSpeciesAddedEvent) {
            ISpecies species = ((SimulationSpeciesAddedEvent)evt).getSpecies();
            for(int i = 0; i < species.getChildTypeCount(); i++) {
                IAtomTypeLeaf newType = species.getChildType(i);
                int indexMax = newType.getIndex();
                agents = Arrays.resizeArray(agents, indexMax+1);
                addAgent(newType);
            }
        }
        else if (evt instanceof SimulationAtomTypeIndexChangedEvent) {
            IAtomTypeLeaf atomType = ((SimulationAtomTypeIndexChangedEvent)evt).getAtomType();
            if (!(atomType instanceof IAtomTypeLeaf)) {
                return;
            }
            int oldIndex = ((SimulationAtomTypeIndexChangedEvent)evt).getOldIndex();
            int newIndex = ((IAtomTypeLeaf)atomType).getIndex();
            if (newIndex >= agents.length) {
                agents = Arrays.resizeArray(agents, newIndex+1);
            }
            agents[newIndex] = agents[oldIndex];
            agents[oldIndex] = null;
        }
        else if (evt instanceof SimulationAtomTypeMaxIndexEvent) {
            int maxIndex = ((SimulationAtomTypeMaxIndexEvent)evt).getMaxIndex();
            agents = Arrays.resizeArray(agents, maxIndex+1);
        }
    }
    
    protected void addAgent(IAtomTypeLeaf type) {
        agents[type.getIndex()] = agentSource.makeAgent(type);
    }
    
    /**
     * Interface for an object that wants an agent associated with each 
     * AtomType in a Simulation.
     */
    public interface AgentSource {
        /**
         * Returns the Class of the agent.  This is used to create an array of 
         * the appropriate Class.
         */
        public Class getSpeciesAgentClass();

        /**
         * Returns an agent for the given AtomType.
         */
        public Object makeAgent(IAtomTypeLeaf type);
        
        /**
         * This informs the agent source that the agent is going away and that 
         * the agent source should disconnect the agent from other elements.
         */
        public void releaseAgent(Object agent, IAtomTypeLeaf type);
    }

    private static final long serialVersionUID = 1L;
    private final AgentSource agentSource;
    protected Object[] agents;
    protected IEventManager simEventManager;
    protected ISpeciesManager speciesManager;
    private final boolean isBackend;

    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator implements Serializable {
        protected AgentIterator(AtomTypeAgentManager agentManager) {
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
        
        private static final long serialVersionUID = 1L;
        private final AtomTypeAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
}
