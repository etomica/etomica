package etomica.atom;

import java.io.Serializable;
import java.lang.reflect.Array;

import etomica.api.IAtomType;
import etomica.simulation.SimulationAtomTypeAddedEvent;
import etomica.simulation.SimulationAtomTypeIndexChangedEvent;
import etomica.simulation.SimulationAtomTypeMaxIndexEvent;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SimulationListener;
import etomica.simulation.SimulationSpeciesRemovedEvent;
import etomica.simulation.SpeciesManager;
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
public class AtomTypeAgentManager implements SimulationListener, java.io.Serializable {

    public AtomTypeAgentManager(AgentSource source) {
        agentSource = source;
        isBackend = true;
    }
    
    public AtomTypeAgentManager(AgentSource source, SpeciesManager speciesManager,
            SimulationEventManager simEventManager, boolean isBackend) {
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
    public void setAgent(AtomType atomType, Object newAgent) {
        agents[atomType.getIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the agent the given AtomType.  For repeated
     * access to the agents from multiple AtomTypes, it might be faster to use 
     * the above getAgents method.
     */
    public Object getAgent(IAtomType type) {
        return agents[type.getIndex()];
    }
    
    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(IAtomType parentType) {
        Object agent = agents[parentType.getIndex()];
        if (agent != null) {
            agentSource.releaseAgent(agent, parentType);
        }
        agents[parentType.getIndex()] = null;

        if (parentType instanceof AtomTypeMolecule) {
            IAtomType[] children = ((AtomTypeMolecule)parentType).getChildTypes();
            for (int i=0; i<children.length; i++) {
                releaseAgents(children[i]);
            }
        }
    }
    
    /**
     * Creates the agents associated with the children of the given AtomType.
     */
    private void makeChildAgents(AtomTypeMolecule parentType) {
        IAtomType[] children = parentType.getChildTypes();
        for (int i=0; i<children.length; i++) {
            addAgent(children[i]);
            if (children[i] instanceof AtomTypeMolecule) {
                makeChildAgents((AtomTypeMolecule)children[i]);
            }
        }
    }

    private void makeAllAgents() {
        AtomTypeMolecule[] moleculeTypes = speciesManager.getMoleculeTypes();
        for (int i=0; i<moleculeTypes.length; i++) {
            addAgent(moleculeTypes[i]);
            makeChildAgents(moleculeTypes[i]);
        }
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private int getGlobalMaxIndex() {
        int max = 0;
        AtomTypeMolecule[] speciesAgentTypes = speciesManager.getMoleculeTypes();
        for (int i=0; i<speciesAgentTypes.length; i++) {
            if (speciesAgentTypes[i].getIndex() > max) {
                max = speciesAgentTypes[i].getIndex();
            }
            int childMax = getMaxIndexOfChildren(speciesAgentTypes[i]);
            if (childMax > max) {
                max = childMax;
            }
        }
        return max;
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private static int getMaxIndexOfChildren(AtomTypeMolecule parentType) {
        int max = 0;
        IAtomType[] children = parentType.getChildTypes();
        for (int i=0; i<children.length; i++) {
            if (children[i].getIndex() > max) {
                max = children[i].getIndex();
            }
            if (children[i] instanceof AtomTypeMolecule) {
                int childMax = getMaxIndexOfChildren((AtomTypeMolecule)children[i]);
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
        // remove ourselves as a listener to the old box
        simEventManager.removeListener(this);
        AtomTypeMolecule[] moleculeTypes = speciesManager.getMoleculeTypes();
        for (int i=0; i<moleculeTypes.length; i++) {
            releaseAgents(moleculeTypes[i]);
        }
        agents = null;
    }
    
    /**
     * Sets the SpeciesRoot for which this AtomAgentManager will manage 
     * AtomType agents.
     */
    public void init(SpeciesManager newSpeciesManager, SimulationEventManager newSimEventManager) {
        simEventManager = newSimEventManager;
        speciesManager = newSpeciesManager;
        simEventManager.addListener(this, isBackend);

        int numTypes = getGlobalMaxIndex()+1;
        
        agents = (Object[])Array.newInstance(agentSource.getTypeAgentClass(), numTypes);
        // fill in the array with agents from all the atoms
        makeAllAgents();
    }
    
    public void actionPerformed(SimulationEvent evt) {
        // we learn about new Species via AtomTypeAdded events
        if (evt instanceof SimulationSpeciesRemovedEvent) {
            AtomTypeMolecule parentType = ((SimulationSpeciesRemovedEvent)evt).getSpecies().getMoleculeType();
            releaseAgents(parentType);
        }
        else if (evt instanceof SimulationAtomTypeAddedEvent) {
            IAtomType newType = ((SimulationAtomTypeAddedEvent)evt).getAtomType();
            int indexMax = newType.getIndex();
            if (newType instanceof AtomTypeMolecule) {
                int childMax = getMaxIndexOfChildren((AtomTypeMolecule)newType);
                if (childMax > indexMax) {
                    indexMax = childMax;
                }
            }
            agents = Arrays.resizeArray(agents, indexMax+1);
            addAgent(newType);
            if (newType instanceof AtomTypeMolecule) {
                makeChildAgents((AtomTypeMolecule)newType);
            }
        }
        else if (evt instanceof SimulationAtomTypeIndexChangedEvent) {
            int oldIndex = ((SimulationAtomTypeIndexChangedEvent)evt).getOldIndex();
            IAtomType atomType = ((SimulationAtomTypeIndexChangedEvent)evt).getAtomType();
            agents[atomType.getIndex()] = agents[oldIndex];
            agents[oldIndex] = null;
        }
        else if (evt instanceof SimulationAtomTypeMaxIndexEvent) {
            int maxIndex = ((SimulationAtomTypeMaxIndexEvent)evt).getMaxIndex();
            agents = Arrays.resizeArray(agents, maxIndex+1);
        }
    }
    
    protected void addAgent(IAtomType type) {
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
        public Class getTypeAgentClass();

        /**
         * Returns an agent for the given AtomType.
         */
        public Object makeAgent(IAtomType type);
        
        /**
         * This informs the agent source that the agent is going away and that 
         * the agent source should disconnect the agent from other elements.
         */
        public void releaseAgent(Object agent, IAtomType type);
    }

    private static final long serialVersionUID = 1L;
    private final AgentSource agentSource;
    protected Object[] agents;
    protected SimulationEventManager simEventManager;
    protected SpeciesManager speciesManager;
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
