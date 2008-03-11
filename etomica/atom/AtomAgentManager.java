package etomica.atom;

import java.io.Serializable;
import java.lang.reflect.Array;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;

import etomica.atom.iterator.AtomIteratorTreeBox;
import etomica.box.BoxAtomAddedEvent;
import etomica.box.BoxAtomEvent;
import etomica.box.BoxAtomIndexChangedEvent;
import etomica.box.BoxAtomRemovedEvent;
import etomica.box.BoxEvent;
import etomica.box.BoxGlobalAtomIndexEvent;
import etomica.box.BoxListener;
import etomica.util.Arrays;

/**
 * AtomAgentManager acts on behalf of client classes (an AgentSource) to manage 
 * agents in every Atom in a box.  When atoms are added or removed from the 
 * box, the agents array (indexed by the atom's global index) is updated.  
 * The client should call getAgents() at any point where an atom might have 
 * have been added to the system because the old array would be stale at that
 * point. 
 * @author andrew
 */
public class AtomAgentManager implements BoxListener, Serializable {

    public AtomAgentManager(AgentSource source, IBox box) {
        this(source, box, true);
    }
    
    public AtomAgentManager(AgentSource source, IBox box, boolean isBackend) {
        agentSource = source;
        this.isBackend = isBackend;
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
        return agents[a.getGlobalIndex()];
    }
    
    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.  The
     * IAtom's old agent is not released.  This should be done manually if
     * needed.
     */
    public void setAgent(IAtom a, Object newAgent) {
        agents[a.getGlobalIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the box the Manager is tracking.
     */
    public IBox getBox(){
        return box;
    }
    
    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        AtomIteratorTreeBox iterator = new AtomIteratorTreeBox(box,Integer.MAX_VALUE,true);
        iterator.reset();
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            // check if atom's spot in the array even exists yet
            if (atom.getGlobalIndex() < agents.length) {
                Object agent = agents[atom.getGlobalIndex()];
                if (agent != null) {
                    agentSource.releaseAgent(agent,atom);
                }
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
                box.getMaxGlobalIndex()+1+reservoirSize);
        // fill in the array with agents from all the atoms
        AtomIteratorTreeBox iterator = new AtomIteratorTreeBox(box,Integer.MAX_VALUE,true);
        iterator.reset();
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            addAgent(atom);
        }
    }
    
    public void actionPerformed(BoxEvent evt) {
        if (evt instanceof BoxAtomEvent) {
            IAtom a = ((BoxAtomEvent)evt).getAtom();
            if (evt instanceof BoxAtomAddedEvent) {
                addAgent(a);
                if (a instanceof IMolecule) {
                    // add all atoms below this atom
                    addAgent(a);
                    IAtomSet childList = ((IMolecule)a).getChildList();
                    for (int i=0; i<childList.getAtomCount(); i++) {
                        addAgent(childList.getAtom(i));
                    }
                }       
            }
            else if (evt instanceof BoxAtomRemovedEvent) {
                int index = a.getGlobalIndex();
                if (agents[index] != null) {
                    // Atom used to have an agent.  nuke it.
                    agentSource.releaseAgent(agents[index], a);
                    agents[index] = null;
                }
                if (a instanceof IMolecule) {
                    // nuke all atoms below this atom
                    IAtomSet childList = ((IMolecule)a).getChildList();
                    for (int i=0; i<childList.getAtomCount(); i++) {
                        IAtom childAtom = childList.getAtom(i);
                        index = childAtom.getGlobalIndex();
                        if (agents[index] != null) {
                            // Atom used to have an agent.  nuke it.
                            agentSource.releaseAgent(agents[index], childAtom);
                            agents[index] = null;
                        }
                    }
                }
            }
            else if (evt instanceof BoxAtomIndexChangedEvent) {
                // the atom's index changed.  assume it would get the same agent
                int oldIndex = ((BoxAtomIndexChangedEvent)evt).getOldIndex();
                agents[a.getGlobalIndex()] = agents[oldIndex];
                agents[oldIndex] = null;
            }
        }
        else if (evt instanceof BoxGlobalAtomIndexEvent) {
            int newMaxIndex = ((BoxGlobalAtomIndexEvent)evt).getMaxIndex();
            if (agents.length > newMaxIndex+reservoirSize || agents.length < newMaxIndex) {
                // indices got compacted.  If our array is a lot bigger than it
                // needs to be, shrink it.
                // ... or we've been notified that atoms are about to get added to the 
                // system.  Make room for them
                agents = Arrays.resizeArray(agents,newMaxIndex+1+reservoirSize);
            }
        }
    }
    
    protected void addAgent(IAtom a) {
        if (agents.length < a.getGlobalIndex()+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = Arrays.resizeArray(agents,a.getGlobalIndex()+1+reservoirSize);
        }
        agents[a.getGlobalIndex()] = agentSource.makeAgent(a);
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
    protected final boolean isBackend;
    protected int reservoirSize;
    
    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator implements Serializable {

        protected AgentIterator(AtomAgentManager agentManager) {
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
        private final AtomAgentManager agentManager;
        private int cursor;
        private Object[] agents;
    }
    
}
