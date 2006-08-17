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
        setRoot(speciesRoot);
        this.isBackend = isBackend;
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
            speciesRoot.getEventManager().removeListener(this);
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
        speciesRoot.getEventManager().addListener(this, isBackend);
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

    private final PhaseAgentSource agentSource;
    protected Object[] agents;
    private SpeciesRoot speciesRoot;
    private final boolean isBackend;
}
