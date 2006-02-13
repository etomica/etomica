package etomica.atom;

import java.lang.reflect.Array;

import etomica.atom.iterator.AtomIteratorTree;
import etomica.phase.Phase;
import etomica.phase.PhaseEvent;
import etomica.phase.PhaseListener;
import etomica.util.Arrays;

/**
 * AtomAgentManager acts on behalf of client classes (an AgentSource) to manage 
 * agents in every Atom in a phase.  When atoms are added or removed from the 
 * phase, the agents array (indexed by the atom's global index) is updated.  
 * The client should call getAgents() at any point where an atom might have 
 * have been added to the system because the old array would be stale at that
 * point. 
 * @author andrew
 */
public class AtomAgentManager implements PhaseListener, java.io.Serializable {

    public AtomAgentManager(AgentSource source, Phase phase) {
        agentSource = source;
        setPhase(phase);
        newAtomList = new AtomList();
    }
    
    public Object[] getAgents() {
        if (!newAtomList.isEmpty()) {
            // if newAtomList has Atoms, agents must have been full, so extend it now
            agents = Arrays.resizeArray(agents,newMaxIndex+phase.getSpeciesMaster().getIndexReservoirSize());
            for (AtomLinker link = newAtomList.header.next; link != newAtomList.header; link = link.next) {
                addAgent(link.atom);
            }
            newAtomList.clear();
        }
        newMaxIndex = -1;
        return agents;
    }
    
    public void setPhase(Phase newPhase) {
        if (phase != null) {
            // remove ourselves as a listener to the old phase
            phase.removeListener(this);
            AtomIteratorTree iterator = new AtomIteratorTree(phase.getSpeciesMaster(),Integer.MAX_VALUE,true);
            iterator.reset();
            while (iterator.hasNext()) {
                Atom atom = iterator.nextAtom();
                agentSource.releaseAgent(agents[atom.getGlobalIndex()],atom);
            }
        }
        phase = newPhase;
        if (phase == null) {
            if (agents != null) {
                // free up the memory
                agents = null;
            }
            return;
        }
        phase.addListener(this);
        // hope the class returns an actual class with a null Atom and use it to construct
        // the array
        agents = (Object[])Array.newInstance(agentSource.getAgentClass(),
                phase.getSpeciesMaster().getMaxGlobalIndex()+1);
        // fill in the array with agents from all the atoms
        AtomIteratorTree iterator = new AtomIteratorTree(phase.getSpeciesMaster(),Integer.MAX_VALUE,true);
        iterator.reset();
        while (iterator.hasNext()) {
            addAgent(iterator.nextAtom());
        }
    }
    
    public void actionPerformed(PhaseEvent evt) {
        if (evt.type() == PhaseEvent.ATOM_ADDED) {
            Atom a = evt.atom();
            if (a.type.isLeaf()) {
                addAgent(a);
            }
            else {
                if (treeIterator == null) {
                    treeIterator = new AtomIteratorTree(Integer.MAX_VALUE);
                    treeIterator.setDoAllNodes(true);
                }
                // add all atoms below this atom
                treeIterator.setRoot(a);
                treeIterator.reset();
                while (treeIterator.hasNext()) {
                    addAgent(treeIterator.nextAtom());
                }
            }       
        }
        else if (evt.type() == PhaseEvent.ATOM_REMOVED) {
            Atom a = evt.atom();
            if (a.type.isLeaf()) {
                int index = a.getGlobalIndex();
                if (agents.length > index && agents[index] != null) {
                    // Atom used to have an agent.  nuke it.
                    agentSource.releaseAgent(agents[index], a);
                    agents[index] = null;
                }
                else if (agents.length < index+1) {
                    // must have been just added
                    newAtomList.remove(a);
                }
            }
            else {
                if (treeIterator == null) {
                    treeIterator = new AtomIteratorTree(Integer.MAX_VALUE);
                    treeIterator.setDoAllNodes(true);
                }
                // nuke all atoms below this atom
                treeIterator.setRoot(a);
                while (treeIterator.hasNext()) {
                    int index = treeIterator.nextAtom().getGlobalIndex();
                    if (agents[index] != null) {
                        // Atom used to have an agent.  nuke it.
                        agents[index] = null;
                    }
                    else if (agents.length < index+1) {
                        // maybe it was just added
                        newAtomList.remove(a);
                    }
                }
            }
        }
        else if (evt.type() == PhaseEvent.ATOM_CHANGE_INDEX) {
            // the atom's index changed.  assume it would get the same agent
            agents[evt.atom().getGlobalIndex()] = agents[evt.getIndex()];
            agents[evt.getIndex()] = null;
        }
        else if (evt.type() == PhaseEvent.GLOBAL_INDEX) {
            // indices got compacted.  If our array is a lot bigger than it
            // needs to be, shrink it.
            SpeciesMaster speciesMaster = (SpeciesMaster)evt.getSource();
            int reservoirSize = speciesMaster.getIndexReservoirSize();
            if (agents.length > evt.getIndex()+reservoirSize) {
                agents = Arrays.resizeArray(agents,evt.getIndex()+reservoirSize);
            }
        }
    }
    
    protected void addAgent(Atom a) {
        if (agents.length < a.getGlobalIndex()+1) {
            // no room in the array.  reallocate the array lazily later.
            newAtomList.add(a);
            if (newMaxIndex < a.getGlobalIndex()) {
                newMaxIndex = a.getGlobalIndex();
            }
            return;
        }
        agents[a.getGlobalIndex()] = agentSource.makeAgent(a);
    }
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface AgentSource {
        public Class getAgentClass();
        
        public Object makeAgent(Atom a);
        
        //allow any agent to be disconnected from other elements 
        public void releaseAgent(Object agent, Atom atom);
    }

    private final AgentSource agentSource;
    protected Object[] agents;
    private AtomList newAtomList;
    private AtomIteratorTree treeIterator;
    private Phase phase;
    private int newMaxIndex;
}
