package etomica.box;

import java.lang.reflect.Array;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SimulationListener;
import etomica.simulation.SimulationBoxAddedEvent;
import etomica.simulation.SimulationBoxEvent;
import etomica.simulation.SimulationBoxRemovedEvent;
import etomica.util.Arrays;

/**
 * BoxAgentManager acts on behalf of client classes (a BoxAgentSource) to manage 
 * agents for each Box in a simulation.  When Box instances are added or removed 
 * from the simulation, the agents array (indexed by the box's index) is updated.  
 * The client should call getAgents() at any point where a Box might have have been 
 * added to (or removed from) the system because the old array would be stale at that
 * point. 
 * @author andrew
 */
public class BoxAgentManager implements SimulationListener, java.io.Serializable {

    public BoxAgentManager(BoxAgentSource source) {
        agentSource = source;
        isBackend = true;
    }

    public BoxAgentManager(BoxAgentSource source, ISimulation sim,
            boolean isBackend) {
        agentSource = source;
        this.isBackend = isBackend;
        setSimulation(sim);
    }
    
    /**
     * Returns the agent associated with the given box
     */
    public Object getAgent(IBox box) {
        return agents[box.getIndex()];
    }
    
    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Sets the Simulation containing Boxs to be tracked.  This method should
     * not be called if setSimulationEventManager is called.
     */
    public void setSimulation(ISimulation sim) {
        simEventManager = sim.getEventManager();
        // this will crash if the given sim is in the middle of its constructor
        simEventManager.addListener(this, isBackend);

        // hope the class returns an actual class with a null Atom and use it to construct
        // the array
        IBox[] boxs = sim.getBoxs();
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),boxs.length);
        for (int i=0; i<boxs.length; i++) {
            addAgent(boxs[i]);
        }
    }
    
    /**
     * Notifies the BoxAgentManager that it should release all agents and 
     * stop listening for events from the simulation.
     */
    public void dispose() {
        // remove ourselves as a listener to the old box
        simEventManager.removeListener(this);
        for (int i=0; i<agents.length; i++) {
            if (agents[i] != null) {
                agentSource.releaseAgent(agents[i]);
            }
        }
        agents = null;
    }
    
    public void actionPerformed(SimulationEvent evt) {
        if (evt instanceof SimulationBoxAddedEvent) {
            addAgent(((SimulationBoxEvent)evt).getBox());
        }
        else if (evt instanceof SimulationBoxRemovedEvent) {
            IBox box = ((SimulationBoxEvent)evt).getBox();
            // The given Box got removed.  The remaining boxs got shifted
            // down.
            int index = box.getIndex();
            agentSource.releaseAgent(agents[index]);
            for (int i=index; i<agents.length-1; i++) {
                agents[i] = agents[i+1];
            }
        }
    }
    
    protected void addAgent(IBox box) {
        agents = Arrays.resizeArray(agents,box.getIndex()+1);
        agents[box.getIndex()] = agentSource.makeAgent(box);
    }
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface BoxAgentSource {
        public Class getAgentClass();
        
        public Object makeAgent(IBox box);
        
        //allow any agent to be disconnected from other elements 
        public void releaseAgent(Object agent); 
    }

    private static final long serialVersionUID = 1L;
    private final BoxAgentSource agentSource;
    protected SimulationEventManager simEventManager;
    protected Object[] agents;
    private final boolean isBackend;
    
    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator {
        protected AgentIterator(BoxAgentManager agentManager) {
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
        
        private final BoxAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
}
