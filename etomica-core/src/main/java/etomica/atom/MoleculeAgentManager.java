/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.box.*;
import etomica.simulation.*;
import etomica.util.Arrays;

import java.io.Serializable;
import java.lang.reflect.Array;

/**
 * MoleculeAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents in every IMolecule in a box.  When molecules are added or
 * removed from the box, the agents array (indexed by the atom's global index)
 * is updated.
 * 
 * @author Andrew Schultz
 */
public class MoleculeAgentManager implements BoxEventListener, SimulationListener, Serializable {

    private static final long serialVersionUID = 1L;
    protected final MoleculeAgentSource agentSource;
    protected final Box box;
    protected final Simulation sim;
    protected Object[][] agents;
    protected int reservoirSize;
    
    public MoleculeAgentManager(Simulation sim, Box box, MoleculeAgentSource source) {
        agentSource = source;
        this.box = box;
        this.sim = sim;
        setReservoirSize(30);
        setupBox();
    }

    public int getReservoirSize() {
        return reservoirSize;
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
    public Object getAgent(IMolecule a) {
        return agents[a.getType().getIndex()][a.getIndex()];
    }
    
    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.  The
     * IAtom's old agent is not released.  This should be done manually if
     * needed.
     */
    public void setAgent(IMolecule a, Object newAgent) {
        agents[a.getType().getIndex()][a.getIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the box the Manager is tracking.
     */
    public Box getBox(){
        return box;
    }
    
    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            IMoleculeList molecules = box.getMoleculeList(sim.getSpecies(i));
            for (int j=0; i<molecules.getMoleculeCount(); j++) {
                // check if atom's spot in the array even exists yet
                IMolecule molecule = molecules.getMolecule(i);
                if (molecule.getIndex() < agents[i].length) {
                    Object agent = agents[i][molecule.getIndex()];
                    if (agent != null) {
                        agentSource.releaseAgent(agent,molecule);
                    }
                }
            }
        }
        agents = null;
    }

    /**
     * Sets the Box in which this AtomAgentManager will manage molecule agents.
     */
    protected void setupBox() {
        sim.getEventManager().addListener(this);
        box.getEventManager().addListener(this);
        agents = new Object[sim.getSpeciesCount()][];
        for (int i=0; i<agents.length; i++) {
            agents[i] = (Object[])Array.newInstance(agentSource.getMoleculeAgentClass(),
                    box.getNMolecules(sim.getSpecies(i)));
        }
        // fill in the array with agents from all the molecules
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            addAgent(molecules.getMolecule(i));
        }
    }
    
    public void boxMoleculeAdded(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        addAgent(mole);
    }
    
    public void boxMoleculeRemoved(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        int index = mole.getIndex();
        int typeIndex = mole.getType().getIndex();
        Object[] speciesAgents = agents[typeIndex];
        if (speciesAgents[index] != null) {
            // Atom used to have an agent.  nuke it.
            agentSource.releaseAgent(speciesAgents[index], mole);
            speciesAgents[index] = null;
        }
    }

    public void boxMoleculeIndexChanged(BoxMoleculeIndexEvent e) {
        IMolecule mole = e.getMolecule();
        // the atom's index changed.  assume it would get the same agent
        int oldIndex = e.getIndex();
        int typeIndex = mole.getType().getIndex();
        Object[] speciesAgents = agents[typeIndex];
        speciesAgents[mole.getIndex()] = speciesAgents[oldIndex];
        speciesAgents[oldIndex] = null;
    }

    public void boxNumberMolecules(BoxMoleculeCountEvent e) {
        int speciesIndex = e.getSpecies().getIndex();
        int newMaxIndex = e.getCount();
        if (agents[speciesIndex].length > newMaxIndex+reservoirSize || agents[speciesIndex].length < newMaxIndex) {
            // indices got compacted.  If our array is a lot bigger than it
            // needs to be, shrink it.
            // ... or we've been notified that atoms are about to get added to the
            // system.  Make room for them
            agents[speciesIndex] = Arrays.resizeArray(agents[speciesIndex],newMaxIndex+reservoirSize);
        }
    }

    public void boxGlobalAtomLeafIndexChanged(BoxIndexEvent e) {}

    public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {}

    public void simulationSpeciesAdded(SimulationSpeciesEvent e) {
        agents = (Object[][])Arrays.resizeArray(agents, agents.length+1);
        agents[agents.length-1] = (Object[])Array.newInstance(agentSource.getMoleculeAgentClass(), 0);
    }

    public void simulationSpeciesRemoved(SimulationSpeciesEvent e) {
        agents = (Object[][])Arrays.removeObject(agents, agents[e.getSpecies().getIndex()]);

    }

    public void simulationBoxAdded(SimulationBoxEvent e) {
    }

    public void simulationBoxRemoved(SimulationBoxEvent e) {
    }

    public void simulationSpeciesIndexChanged(SimulationSpeciesIndexEvent e) {
    }

    public void simulationSpeciesMaxIndexChanged(SimulationIndexEvent e) {
    }

    public void simulationAtomTypeIndexChanged(SimulationAtomTypeEvent e) {
    }

    public void simulationAtomTypeMaxIndexChanged(SimulationIndexEvent e) {
    }

    protected void addAgent(IMolecule a) {
        int speciesIndex = a.getType().getIndex();
        if (agents[speciesIndex].length < a.getIndex()+1) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents[speciesIndex] = Arrays.resizeArray(agents[speciesIndex],a.getIndex()+1+reservoirSize);
        }
        agents[speciesIndex][a.getIndex()] = agentSource.makeAgent(a);
    }
    /**
     * Interface for an object that wants an agent associated with each
     * IMolecule a Box.
     */
    public interface MoleculeAgentSource {
        /**
         * Returns the Class of the agent.  This is used to create an array of
         * the appropriate Class.
         */
        Class getMoleculeAgentClass();

        /**
         * Returns an agent for the given Atom.
         */
        Object makeAgent(IMolecule a);

        /**
         * This informs the agent source that the agent is going away and that
         * the agent source should disconnect the agent from other elements
         */
        void releaseAgent(Object agent, IMolecule atom);
    }
    
    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator implements Serializable {

        private static final long serialVersionUID = 1L;
        private final MoleculeAgentManager agentManager;
        private int speciesCursor, moleculeCursor;
        private Object[][] agents;
        private Object[] speciesAgents;
        protected AgentIterator(MoleculeAgentManager agentManager) {
            this.agentManager = agentManager;
        }

        public void reset() {
            speciesCursor = 0;
            moleculeCursor = 0;
            agents = agentManager.agents;
            if (agents.length > 0) {
                speciesAgents = agents[0];
            }
            else {
                speciesAgents = null;
            }
        }

        public boolean hasNext() {
            if (speciesAgents == null) {
                return false;
            }
            do {
                while (moleculeCursor == speciesAgents.length) {
                    moleculeCursor = 0;
                    speciesCursor++;
                    if (speciesCursor < agents.length) {
                        speciesAgents = agents[speciesCursor];
                    }
                    else {
                        speciesAgents = null;
                        return false;
                    }
                }
                while (moleculeCursor < speciesAgents.length) {
                    if (speciesAgents[moleculeCursor] != null) {
                        return true;
                    }
                    moleculeCursor++;
                }
            } while (true);
        }

        public Object next() {
            if (speciesAgents == null) {
                return false;
            }
            do {
                while (moleculeCursor == speciesAgents.length) {
                    moleculeCursor = 0;
                    if (speciesCursor < agents.length) {
                        speciesCursor++;
                        speciesAgents = agents[speciesCursor];
                    }
                    else {
                        speciesAgents = null;
                        return null;
                    }
                }
                while (moleculeCursor < speciesAgents.length) {
                    if (speciesAgents[moleculeCursor] != null) {
                        return speciesAgents[moleculeCursor++];
                    }
                }
            } while (true);
        }
    }
}
