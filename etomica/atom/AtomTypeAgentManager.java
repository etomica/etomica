package etomica.atom;

import java.lang.reflect.Array;

import etomica.simulation.SimulationAtomTypeAddedEvent;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;
import etomica.simulation.SimulationSpeciesAddedEvent;
import etomica.simulation.SimulationSpeciesEvent;
import etomica.simulation.SimulationSpeciesRemovedEvent;
import etomica.util.Arrays;

/**
 * AtomTypeAgentManager acts on behalf of client classes (an AgentSource) to 
 * manage agents in every AtomType in a phase.  When species are added or 
 * removed from the simulation, the agents array (indexed by the AtomType's 
 * global index) is updated.  The client should call getAgents() at any point 
 * where an atom might have have been added to the system because the old array
 * would be stale at that point.
 * @author andrew
 */
public class AtomTypeAgentManager implements SimulationListener, java.io.Serializable {

    public AtomTypeAgentManager(AgentSource source, SpeciesRoot root) {
        this(source, root, true);
    }
    
    public AtomTypeAgentManager(AgentSource source, SpeciesRoot root, boolean isBackend) {
        agentSource = source;
        this.isBackend = isBackend;
        if (root != null) {
            setRoot(root);
        }
    }        
    
    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Returns the array of AtomType agents, indexed by the AtomType's index.  
     * The array is of the type returned by the AgentSource's getAgentClass 
     * method.
     */
    public Object[] getAgents() {
        return agents;
    }
    
    /**
     * Convenience method to return the agent the given AtomType.  For repeated
     * access to the agents from multiple AtomTypes, it might be faster to use 
     * the above getAgents method.
     */
    public Object getAgent(AtomType type) {
        return agents[type.getIndex()];
    }
    
    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(AtomType parentType) {
        Object agent = agents[parentType.getIndex()];
        if (agent != null) {
            agentSource.releaseAgent(agent, parentType);
        }
        agents[parentType.getIndex()] = null;

        if (parentType instanceof AtomTypeGroup) {
            AtomType[] children = ((AtomTypeGroup)parentType).getChildTypes();
            for (int i=0; i<children.length; i++) {
                if (children[i] instanceof AtomTypeGroup) {
                    releaseAgents(children[i]);
                }
            }
        }
    }
    
    /**
     * Creates the agents associated with the children of the given AtomType.
     */
    private void makeChildAgents(AtomTypeGroup parentType) {
        AtomType[] children = parentType.getChildTypes();
        for (int i=0; i<children.length; i++) {
            addAgent(children[i]);
            if (children[i] instanceof AtomTypeGroup) {
                makeChildAgents((AtomTypeGroup)children[i]);
            }
        }
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private static int getMaxIndexOfChildren(AtomTypeGroup parentType) {
        int max = 0;
        AtomType[] children = parentType.getChildTypes();
        for (int i=0; i<children.length; i++) {
            if (children[i].getIndex() > max) {
                max = children[i].getIndex();
            }
            if (children[i] instanceof AtomTypeGroup) {
                int childMax = getMaxIndexOfChildren((AtomTypeGroup)children[i]);
                if (childMax > max) {
                    max = childMax;
                }
            }
        }
        return max;
    }
    
    /**
     * Unregisters this class as a listener for AtomType-related events and 
     * releases its agents.
     */
    public void dispose() {
        // remove ourselves as a listener to the old phase
        root.getEventManager().removeListener(this);
        releaseAgents(root.type);
        agents = null;
    }
    
    /**
     * Sets the SpeciesRoot for which this AtomAgentManager will manage 
     * AtomType agents.
     */
    public void setRoot(SpeciesRoot newRoot) {
        root = newRoot;
        root.getEventManager().addListener(this, isBackend);

        int numTypes = getMaxIndexOfChildren((AtomTypeGroup)root.type)+1;
        
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),
                numTypes);
        // fill in the array with agents from all the atoms
        addAgent(root.type);
        makeChildAgents((AtomTypeGroup)root.type);
    }
    
    public void actionPerformed(SimulationEvent evt) {
        // we learn about new Species via AtomTypeAdded events
        if (evt instanceof SimulationSpeciesRemovedEvent) {
            AtomTypeGroup parentType = ((SimulationSpeciesRemovedEvent)evt).getSpecies().getMoleculeType().getParentType();
            releaseAgents(parentType);
        }
        else if (evt instanceof SimulationAtomTypeAddedEvent) {
            AtomType newType = ((SimulationAtomTypeAddedEvent)evt).getAtomType();
            AtomTypeGroup parentType = newType.getParentType();
            int childMax = getMaxIndexOfChildren(parentType);
            agents = Arrays.resizeArray(agents, childMax+1);
            addAgent(newType);
            if (newType instanceof AtomTypeGroup) {
                makeChildAgents((AtomTypeGroup)newType);
            }
        }
        else if (evt instanceof SimulationAtomTypeCompactedEvent) {
            int typeStart = ((SimulationAtomTypeCompactedEvent)evt).getStartIndex(); // first AtomType removed
            int typeStop = ((SimulationAtomTypeCompactedEvent)evt).getStopIndex();   // first AtomType not removed
            // shift the agents back
            if (typeStop > -1) {
                // typeStop=-1 means the last AtomType removed was the last AtomType
                for (int i=0; i<agents.length-typeStop && i < typeStop-typeStart; i++) {
                    agents[typeStart+i] = agents[typeStop+i];
                }
            }
            // resize the array to drop the extra agent elements
            agents = Arrays.resizeArray(agents, agents.length-(typeStop-typeStart));
        }
    }
    
    protected void addAgent(AtomType type) {
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
        public Class getAgentClass();

        /**
         * Returns an agent for the given AtomType.
         */
        public Object makeAgent(AtomType type);
        
        /**
         * This informs the agent source that the agent is going away and that 
         * the agent source should disconnect the agent from other elements.
         */
        public void releaseAgent(Object agent, AtomType type);
    }

    private static final long serialVersionUID = 1L;
    private final AgentSource agentSource;
    protected Object[] agents;
    private SpeciesRoot root;
    private final boolean isBackend;

    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator {
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
        
        private final AtomTypeAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
}
