package etomica.atom;

import java.lang.reflect.Array;

import etomica.box.Box;
import etomica.box.BoxAtomAddedEvent;
import etomica.box.BoxAtomEvent;
import etomica.box.BoxAtomLeafIndexChangedEvent;
import etomica.box.BoxAtomRemovedEvent;
import etomica.box.BoxEvent;
import etomica.box.BoxGlobalAtomLeafIndexEvent;
import etomica.util.Arrays;

/**
 * AtomLeafAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents for every leaf Atom in a box.  When leaf atoms are added or
 * removed from the box, the agents array (indexed by the atom's global
 * index) is updated.  The client can access and modify the agents via getAgent
 * and setAgent.
 * @author Andrew Schultz
 */
public class AtomLeafAgentManager extends AtomAgentManager {

    public AtomLeafAgentManager(AgentSource source, Box box) {
        this(source, box, true);
    }
    
    public AtomLeafAgentManager(AgentSource source, Box box, boolean isBackend) {
        super(source, box, isBackend);
        // we just want the leaf atoms
    }        
    
    /**
     * Returns the agent associated with the given IAtom.  The IAtom must be
     * from the Box associated with this instance.
     */
    public Object getAgent(IAtom a) {
        return agents[box.getLeafIndex(a)];
    }
    
    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.
     */
    public void setAgent(IAtom a, Object newAgent) {
        agents[box.getLeafIndex(a)] = newAgent;
    }
    
    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        AtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            // leaf index corresponds to the position in the leaf list
            Object agent = agents[i];
            if (agent != null) {
                agentSource.releaseAgent(agent,leafList.getAtom(i));
            }
        }
        agents = null;
    }
    
    /**
     * Sets the Box in which this AtomAgentManager will manage Atom agents.
     */
    protected void setupBox() {
        box.getEventManager().addListener(this, isBackend);
        
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),
                box.getLeafList().getAtomCount()+1+box.getIndexReservoirSize());
        // fill in the array with agents from all the atoms
        AtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            // leaf list position is the leaf index, so don't bother looking
            // that up again.
           addAgent(leafList.getAtom(i), i);
        }
    }
    
    public void actionPerformed(BoxEvent evt) {
        if (evt instanceof BoxAtomEvent) {
            IAtom a = ((BoxAtomEvent)evt).getAtom();
            if (evt instanceof BoxAtomAddedEvent) {
                if (a instanceof IMolecule) {
                    // add all leaf atoms below this atom
                    AtomSet childList = ((IMolecule)a).getChildList();
                    for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
                        addAgent(childList.getAtom(iChild));
                    }
                }
                else {
                    // the atom itself is a leaf
                    addAgent(a);
                }
            }
            else if (evt instanceof BoxAtomRemovedEvent) {
                if (a instanceof IMolecule) {
                    // IAtomGroups don't have agents, but nuke all atoms below this atom
                    AtomSet childList = ((IMolecule)a).getChildList();
                    for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
                        IAtom childAtom = childList.getAtom(iChild);
                        int index = box.getLeafIndex(childAtom);
                        if (agents[index] != null) {
                            // Atom used to have an agent.  nuke it.
                            agentSource.releaseAgent(agents[index], childAtom);
                            agents[index] = null;
                        }
                    }
                }
                else {
                    int index = box.getLeafIndex(a);
                    if (agents[index] != null) {
                        // Atom used to have an agent.  nuke it.
                        agentSource.releaseAgent(agents[index], a);
                        agents[index] = null;
                    }
                }
            }
            else if (evt instanceof BoxAtomLeafIndexChangedEvent) {
                // the atom's index changed.  assume it would get the same agent
                int oldIndex = ((BoxAtomLeafIndexChangedEvent)evt).getOldIndex();
                agents[box.getLeafIndex(a)] = agents[oldIndex];
                agents[oldIndex] = null;
            }
        }
        else if (evt instanceof BoxGlobalAtomLeafIndexEvent) {
            int reservoirSize = box.getIndexReservoirSize();
            // don't use leafList.size() since the SpeciesMaster might be notifying
            // us that it's about to add leaf atoms
            int newMaxIndex = ((BoxGlobalAtomLeafIndexEvent)evt).getMaxIndex();
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
        addAgent(a, box.getLeafIndex(a));
    }
    
    /**
     * Adds an agent for the given leaf atom to the agents array at the given
     * index.
     */
    protected void addAgent(IAtom a, int index) {
        if (agents.length < index+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,index+1+box.getIndexReservoirSize());
        }
        agents[index] = agentSource.makeAgent(a);
    }        
    
    private static final long serialVersionUID = 1L;
    
}
