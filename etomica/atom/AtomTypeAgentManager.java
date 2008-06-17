package etomica.atom;

import java.io.Serializable;
import java.lang.reflect.Array;

import etomica.api.IAtomType;
import etomica.api.ISimulationEventManager;
import etomica.api.ISpecies;
import etomica.api.ISpeciesManager;
import etomica.simulation.SimulationAtomTypeAddedEvent;
import etomica.simulation.SimulationAtomTypeIndexChangedEvent;
import etomica.simulation.SimulationAtomTypeMaxIndexEvent;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;
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
public class AtomTypeAgentManager implements SimulationListener, java.io.Serializable {

    public AtomTypeAgentManager(AgentSource source) {
        agentSource = source;
        isBackend = true;
    }
    
    public AtomTypeAgentManager(AgentSource source, ISpeciesManager speciesManager,
            ISimulationEventManager simEventManager, boolean isBackend) {
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
    public void setAgent(IAtomType atomType, Object newAgent) {
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

        if (parentType instanceof ISpecies) {
            for (int i=0; i<((ISpecies)parentType).getChildTypeCount(); i++) {
                releaseAgents(((ISpecies)parentType).getChildType(i));
            }
        }
    }
    
    /**
     * Creates the agents associated with the children of the given AtomType.
     */
    private void makeChildAgents(ISpecies parentType) {
        for (int i=0; i<parentType.getChildTypeCount(); i++) {
            addAgent(parentType.getChildType(i));
            if (parentType.getChildType(i) instanceof ISpecies) {
                makeChildAgents((ISpecies)parentType.getChildType(i));
            }
        }
    }

    private void makeAllAgents() {
        ISpecies[] moleculeTypes = speciesManager.getMoleculeTypes();
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
        ISpecies[] speciesAgentTypes = speciesManager.getMoleculeTypes();
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
    private static int getMaxIndexOfChildren(ISpecies parentType) {
        int max = 0;
        for (int i=0; i<parentType.getChildTypeCount(); i++) {
            if (parentType.getChildType(i).getIndex() > max) {
                max = parentType.getChildType(i).getIndex();
            }
            if (parentType.getChildType(i) instanceof ISpecies) {
                int childMax = getMaxIndexOfChildren((ISpecies)parentType.getChildType(i));
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
        ISpecies[] moleculeTypes = speciesManager.getMoleculeTypes();
        for (int i=0; i<moleculeTypes.length; i++) {
            releaseAgents(moleculeTypes[i]);
        }
        agents = null;
    }
    
    /**
     * Sets the SpeciesRoot for which this AtomAgentManager will manage 
     * AtomType agents.
     */
    public void init(ISpeciesManager newSpeciesManager, ISimulationEventManager newSimEventManager) {
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
            releaseAgents(((SimulationSpeciesRemovedEvent)evt).getSpecies());
        }
        else if (evt instanceof SimulationAtomTypeAddedEvent) {
            IAtomType newType = ((SimulationAtomTypeAddedEvent)evt).getAtomType();
            int indexMax = newType.getIndex();
            if (newType instanceof ISpecies) {
                int childMax = getMaxIndexOfChildren((ISpecies)newType);
                if (childMax > indexMax) {
                    indexMax = childMax;
                }
            }
            agents = Arrays.resizeArray(agents, indexMax+1);
            addAgent(newType);
            if (newType instanceof ISpecies) {
                makeChildAgents((ISpecies)newType);
            }
        }
        else if (evt instanceof SimulationAtomTypeIndexChangedEvent) {
            int oldIndex = ((SimulationAtomTypeIndexChangedEvent)evt).getOldIndex();
            IAtomType atomType = ((SimulationAtomTypeIndexChangedEvent)evt).getAtomType();
            if (atomType.getIndex() >= agents.length) {
                agents = Arrays.resizeArray(agents, atomType.getIndex()+1);
            }
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
    protected ISimulationEventManager simEventManager;
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
