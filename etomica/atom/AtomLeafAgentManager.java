package etomica.atom;

import java.io.Serializable;
import java.lang.reflect.Array;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IBoxAtomIndexEvent;
import etomica.api.IBoxIndexEvent;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;
import etomica.box.BoxListenerAdapter;
import etomica.util.Arrays;

/**
 * AtomAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents for every leaf Atom in a box.  When leaf atoms are added or
 * removed from the box, the agents array (indexed by the atom's global
 * index) is updated.  The client can access and modify the agents via getAgent
 * and setAgent.
 * 
 * @author Andrew Schultz
 */
public class AtomLeafAgentManager extends BoxListenerAdapter implements Serializable {
    
    public AtomLeafAgentManager(AgentSource source, IBox box) {
        agentSource = source;
        this.box = box;
        setReservoirSize(30);
        setupBox();
    }
    
    /**
     * Sets the size of the manager's "reservoir".  When an atom is removed,
     * the agents array will only be trimmed if the number of holes in the
     * array exceeds the reservoir size.  Also, when the array has no holes and
     * another atom is added, the array will resized to be 
     * numAtoms+reservoirSize to avoid reallocating a new array every time an
     * atom is added.  reservoirSize=0 means the array will
     * always be the same size as the number of atoms (no holes).
     * 
     * The default reservoir size is 30.
     */
    public void setReservoirSize(int newReservoirSize) {
        reservoirSize = newReservoirSize;
    }

    public int getReservoirSize() {
        return reservoirSize;
    }

    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Returns the agent associated with the given IAtom.  The IAtom must be
     * from the Box associated with this instance.
     */
    public Object getAgent(IAtom a) {
        int idx = a.getLeafIndex();
        if (idx < agents.length) {
            return agents[idx];
        }
        return null;
    }
    
    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.  The
     * IAtom's old agent is not released.  This should be done manually if
     * needed.
     */
    public void setAgent(IAtom a, Object newAgent) {
        int idx = a.getLeafIndex();
        if (idx >= agents.length) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,idx+1+reservoirSize);
        }
        agents[a.getLeafIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the box the Manager is tracking.
     */
    public IBox getBox(){
        return box;
    }
    
    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener and release all of its agents.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            // leaf index corresponds to the position in the leaf list
            Object agent = agents[i];
            if (agent != null && agentSource != null) {
                agentSource.releaseAgent(agent, leafList.getAtom(i));
            }
        }
        agents = null;
    }
    
    /**
     * Sets the Box in which this AtomAgentManager will manage Atom agents.
     */
    protected void setupBox() {
        box.getEventManager().addListener(this);
        
        agents = (Object[])Array.newInstance(agentSource != null ? agentSource.getAgentClass() : Object.class,
                box.getLeafList().getAtomCount()+1+reservoirSize);
        // fill in the array with agents from all the atoms
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            // leaf list position is the leaf index, so don't bother looking
            // that up again.
           addAgent(leafList.getAtom(i), i);
        }
    }
    
    public void boxMoleculeAdded(IBoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        // add all leaf atoms below this atom
        IAtomList childList = mole.getChildList();

        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            addAgent(childList.getAtom(iChild));
        }
    }
    
    public void boxMoleculeRemoved(IBoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        // IAtomGroups don't have agents, but nuke all atoms below this atom
        IAtomList childList = mole.getChildList();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtom childAtom = childList.getAtom(iChild);
            int index = childAtom.getLeafIndex();
            if (agents[index] != null) {
                // Atom used to have an agent.  nuke it.
                if (agentSource != null) {
                    agentSource.releaseAgent(agents[index], childAtom);
                }
                agents[index] = null;
            }
        }
    }
    
    public void boxAtomLeafIndexChanged(IBoxAtomIndexEvent e) {
        IAtom a = e.getAtom();
        // the atom's index changed.  assume it would get the same agent
        int oldIndex = e.getIndex();
        agents[((IAtom)a).getLeafIndex()] = agents[oldIndex];
        agents[oldIndex] = null;
    }
    
    public void boxGlobalAtomLeafIndexChanged(IBoxIndexEvent e) {
        // don't use leafList.size() since the SpeciesMaster might be notifying
        // us that it's about to add leaf atoms
        int newMaxIndex = e.getIndex();
        if (agents.length > newMaxIndex+reservoirSize || agents.length < newMaxIndex) {
            // indices got compacted.  If our array is a lot bigger than it
            // needs to be, shrink it.
            // ... or we've been notified that atoms are about to get added to the 
            // system.  Make room for them
            agents = Arrays.resizeArray(agents,newMaxIndex+1+reservoirSize);
        }
    }
    
    /**
     * Adds an agent for the given leaf atom to the agents array.
     */
    protected void addAgent(IAtom a) {
        addAgent(a, a.getLeafIndex());
    }
    
    /**
     * Adds an agent for the given leaf atom to the agents array at the given
     * index.
     */
    protected void addAgent(IAtom a, int index) {
        if (agents.length < index+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,index+1+reservoirSize);
        }
        if (agentSource != null) {
            agents[index] = agentSource.makeAgent(a);
        }
    }        
    
    /**
     * Interface for an object that wants an agent associated with each Atom in
     * a Box.
     */
    public interface AgentSource {
        /**
         * Returns the Class of the agent.  This is used to create an array of 
         * the appropriate Class.
         */
        public Class getAgentClass();

        /**
         * Returns an agent for the given Atom.
         */
        public Object makeAgent(IAtom a);
        
        /**
         * This informs the agent source that the agent is going away and that 
         * the agent source should disconnect the agent from other elements
         */
        public void releaseAgent(Object agent, IAtom atom);
    }

    private static final long serialVersionUID = 1L;
    protected final AgentSource agentSource;
    protected Object[] agents;
    protected final IBox box;
    protected int reservoirSize;
    
    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator implements Serializable {

        protected AgentIterator(AtomLeafAgentManager agentManager) {
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
        private final AtomLeafAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
    
}
