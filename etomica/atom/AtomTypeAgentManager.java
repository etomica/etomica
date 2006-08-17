package etomica.atom;

import java.lang.reflect.Array;

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
        setRoot(root);
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
     * Creates the agents associated with the given AtomType and its children.
     */
    private void makeAgents(AtomType parentType) {
        addAgent(parentType);
        
        if (parentType instanceof AtomTypeGroup) {
            AtomType[] children = ((AtomTypeGroup)parentType).getChildTypes();
            for (int i=0; i<children.length; i++) {
                addAgent(children[i]);
                if (children[i] instanceof AtomTypeGroup) {
                    makeAgents(children[i]);
                }
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
     * Sets the SpeciesRoot for which this AtomAgentManager will manage 
     * AtomType agents.  Setting the root to null will notify the 
     * AtomAgentManager it should disconnect itself as a listener and release
     * all agents.
     */
    public void setRoot(SpeciesRoot newRoot) {
        if (root!= null) {
            // remove ourselves as a listener to the old phase
            root.getEventManager().removeListener(this);
            releaseAgents(root.type);
        }
        root = newRoot;
        if (root == null) {
            agents = null;
            return;
        }
        root.getEventManager().addListener(this, isBackend);

        int numTypes = getMaxIndexOfChildren((AtomTypeGroup)root.type);
        
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),
                numTypes);
        // fill in the array with agents from all the atoms
        makeAgents(root.type);
    }
    
    public void actionPerformed(SimulationEvent evt) {
        if (evt instanceof SimulationSpeciesAddedEvent) {
            AtomTypeGroup parentType = ((SimulationSpeciesEvent)evt).getSpecies().getMoleculeType().getParentType();
            int childMax = getMaxIndexOfChildren(parentType);
            agents = Arrays.resizeArray(agents, childMax);
            makeAgents(parentType);
        }
        else if (evt instanceof SimulationSpeciesRemovedEvent) {
            AtomTypeGroup parentType = ((SimulationSpeciesRemovedEvent)evt).getSpecies().getMoleculeType().getParentType();
            releaseAgents(parentType);
        }
        // handle other events for AtomType index changing and compaction
    }
    
    protected void addAgent(AtomType type) {
        if (agents.length < type.getIndex()+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,type.getIndex()+1);
        }
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
}
