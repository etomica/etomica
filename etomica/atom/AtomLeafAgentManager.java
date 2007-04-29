package etomica.atom;

import java.lang.reflect.Array;

import etomica.phase.Phase;
import etomica.phase.PhaseAtomAddedEvent;
import etomica.phase.PhaseAtomEvent;
import etomica.phase.PhaseAtomLeafIndexChangedEvent;
import etomica.phase.PhaseAtomRemovedEvent;
import etomica.phase.PhaseEvent;
import etomica.phase.PhaseGlobalAtomLeafIndexEvent;
import etomica.util.Arrays;

/**
 * AtomLeafAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents for every leaf Atom in a phase.  When leaf atoms are added or
 * removed from the phase, the agents array (indexed by the atom's global
 * index) is updated.  The client can access and modify the agents via getAgent
 * and setAgent.
 * @author Andrew Schultz
 */
public class AtomLeafAgentManager extends AtomAgentManager {

    public AtomLeafAgentManager(AgentSource source, Phase phase) {
        this(source, phase, true);
    }
    
    public AtomLeafAgentManager(AgentSource source, Phase phase, boolean isBackend) {
        super(source, phase, isBackend);
        // we just want the leaf atoms
        treeIterator.setDoAllNodes(false);
    }        
    
    /**
     * Returns the agent associated with the given IAtom.  The IAtom must be
     * from the Phase associated with this instance.
     */
    public Object getAgent(IAtom a) {
        return agents[speciesMaster.getLeafIndex(a)];
    }
    
    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Phase associated with this instance.
     */
    public void setAgent(IAtom a, Object newAgent) {
        agents[speciesMaster.getLeafIndex(a)] = newAgent;
    }
    
    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener.
     */
    public void dispose() {
        // remove ourselves as a listener to the phase
        speciesMaster.getPhase().getEventManager().removeListener(this);
        AtomArrayList leafList = speciesMaster.getLeafList();
        int nLeaf = leafList.size();
        for (int i=0; i<nLeaf; i++) {
            // leaf index corresponds to the position in the leaf list
            Object agent = agents[i];
            if (agent != null) {
                agentSource.releaseAgent(agent,leafList.get(i));
            }
        }
        agents = null;
    }
    
    /**
     * Sets the Phase in which this AtomAgentManager will manage Atom agents.
     */
    protected void setupPhase() {
        speciesMaster.getPhase().getEventManager().addListener(this, isBackend);
        
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),
                speciesMaster.getLeafList().size()+1+speciesMaster.getIndexReservoirSize());
        // fill in the array with agents from all the atoms
        AtomArrayList leafList = speciesMaster.getLeafList();
        int nLeaf = leafList.size();
        for (int i=0; i<nLeaf; i++) {
            // leaf list position is the leaf index, so don't bother looking
            // that up again.
           addAgent(leafList.get(i), i);
        }
    }
    
    public void actionPerformed(PhaseEvent evt) {
        if (evt instanceof PhaseAtomEvent) {
            IAtom a = ((PhaseAtomEvent)evt).getAtom();
            if (evt instanceof PhaseAtomAddedEvent) {
                if (a instanceof IAtomGroup) {
                    // add all leaf atoms below this atom
                    treeIterator.setRootAtom(a);
                    treeIterator.reset();
                    
                    for (IAtom atom = treeIterator.nextAtom(); atom != null; atom = treeIterator.nextAtom()) {
                        addAgent(atom);
                    }
                }
                else {
                    // the atom itself is a leaf
                    addAgent(a);
                }
            }
            else if (evt instanceof PhaseAtomRemovedEvent) {
                if (a instanceof IAtomGroup) {
                    // IAtomGroups don't have agents, but nuke all atoms below this atom
                    treeIterator.setRootAtom(a);
                    treeIterator.reset();
                    for (IAtom childAtom = treeIterator.nextAtom(); childAtom != null; childAtom = treeIterator.nextAtom()) {
                        int index = speciesMaster.getLeafIndex(childAtom);
                        if (agents[index] != null) {
                            // Atom used to have an agent.  nuke it.
                            agentSource.releaseAgent(agents[index], childAtom);
                            agents[index] = null;
                        }
                    }
                }
                else {
                    int index = speciesMaster.getLeafIndex(a);
                    if (agents[index] != null) {
                        // Atom used to have an agent.  nuke it.
                        agentSource.releaseAgent(agents[index], a);
                        agents[index] = null;
                    }
                }
            }
            else if (evt instanceof PhaseAtomLeafIndexChangedEvent) {
                // the atom's index changed.  assume it would get the same agent
                int oldIndex = ((PhaseAtomLeafIndexChangedEvent)evt).getOldIndex();
                agents[speciesMaster.getLeafIndex(a)] = agents[oldIndex];
                agents[oldIndex] = null;
            }
        }
        else if (evt instanceof PhaseGlobalAtomLeafIndexEvent) {
            int reservoirSize = speciesMaster.getIndexReservoirSize();
            // don't use leafList.size() since the SpeciesMaster might be notifying
            // us that it's about to add leaf atoms
            int newMaxIndex = ((PhaseGlobalAtomLeafIndexEvent)evt).getMaxIndex();
            if (agents.length > newMaxIndex+reservoirSize || agents.length < newMaxIndex) {
                // indices got compacted.  If our array is a lot bigger than it
                // needs to be, shrink it.
                // ... or we've been notified that atoms are about to get added to the 
                // system.  Make room for them
                agents = Arrays.resizeArray(agents,newMaxIndex+1+reservoirSize);
            }
        }
    }

    /**
     * Adds an agent for the given leaf atom to the agents array.
     */
    protected void addAgent(IAtom a) {
        addAgent(a, speciesMaster.getLeafIndex(a));
    }
    
    /**
     * Adds an agent for the given leaf atom to the agents array at the given
     * index.
     */
    protected void addAgent(IAtom a, int index) {
        if (agents.length < index+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,index+1+speciesMaster.getIndexReservoirSize());
        }
        agents[index] = agentSource.makeAgent(a);
    }        
    
    private static final long serialVersionUID = 1L;
    
}
