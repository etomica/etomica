package etomica.atom;

import etomica.atom.iterator.AtomIteratorTreePhase;
import etomica.atom.iterator.AtomIteratorTreeRoot;
import etomica.phase.Phase;
import etomica.phase.PhaseAtomAddedEvent;
import etomica.phase.PhaseAtomIndexChangedEvent;
import etomica.phase.PhaseAtomRemovedEvent;
import etomica.phase.PhaseEventManager;
import etomica.phase.PhaseGlobalAtomIndexEvent;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

/**
 * Coordinator of all species agents in a phase. Parent is SpeciesRoot, and
 * children are SpeciesAgent instances. All instances of SpeciesMaster in a
 * given simulation share the same AtomType.
 * 
 * @author David Kofke and Andrew Schultz
 */

public final class SpeciesMaster implements java.io.Serializable {

    public SpeciesMaster(Phase p, PhaseEventManager eventManager) {
        agentList = new AtomArrayList();
        phaseEventManager = eventManager;
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        reservoirCount = 0;
        phase = p;
        treeIteratorRoot = new AtomIteratorTreeRoot();
        treeIteratorRoot.setDoAllNodes(true);
        treeIteratorRoot.setIterationDepth(Integer.MAX_VALUE);
        treeIteratorPhase = new AtomIteratorTreePhase();
    }

    /**
     * Adds the given SpeciesAgent to the SpeciesMaster's list of
     * SpeciesAgents.  This method should be called by Species.
     */
    public void addSpeciesAgent(SpeciesAgent newSpeciesAgent) {
        newSpeciesAgent.setIndex(agentList.size());
        agentList.add(newSpeciesAgent);
        addAtomNotify(newSpeciesAgent);
    }
    
    /**
     * Notifies the SpeciesMaster that a Species has been removed.  This method
     * should only be called by the SpeciesManager.
     */
    public void removeSpeciesNotify(Species species) {
        for (int i=0; i<agentList.size(); i++) {
            IAtom speciesAgent = agentList.get(i);
            if (speciesAgent.getType().getSpecies() == species) {
                agentList.removeAndReplace(i);
                if (agentList.size() > i) {
                    agentList.get(i).setIndex(i);
                }
                agentList.maybeTrimToSize();
                removeAtomNotify(speciesAgent);
            }
        }
    }

    public AtomArrayList getAgentList() {
        return agentList;
    }

    /**
     * Returns the number of molecules in the Phase
     */
    public int moleculeCount() {
        return moleculeCount;
    }

    /**
     * Returns an AtomArrayList containing the leaf atoms in the Phase
     */
    public AtomArrayList getLeafList() {
        return leafList;
    }
    
    public Phase getPhase() {
        return phase;
    }

    /**
     * Returns a "global" index for the Phase.  This method should only be
     * called by Atom.
     */
    public int requestGlobalIndex() {
        if (reservoirCount == 0) {
            return ++maxIndex;
        }
        reservoirCount--;
        return indexReservoir[reservoirCount];
    }
    
    protected void returnGlobalIndex(int index) {
        // add the index to the reservoir first
        indexReservoir[reservoirCount] = index;
        reservoirCount++;
        if (reservoirCount == reservoirSize) {
            // reservoir is full
            collapseGlobalIndices();
        }
    }
    
    /**
     * Returns the maximum global index for the Phase.
     */
    public int getMaxGlobalIndex() {
        return maxIndex;
    }
    
    /**
     * Collapses the global indices.  At the end, the maximum global index of 
     * an atom in the simulation will be maxIndex.  At the start, if indices 
     * in the reservoir are greater than the final maximum global index, they
     * are discarded.  Then Atom with global indices greater than the final 
     * maximum global index are given indices from the reservoir that are less
     * than the final maximum global index.  At the end maxIndex is decreased 
     * by maxIndex.
     */
    private void collapseGlobalIndices() {
        for (int j=0; j<reservoirCount; ) {
            if (indexReservoir[j] > maxIndex-reservoirSize) {
                // this index isn't useful to us, so just drop it
                reservoirCount--;
                indexReservoir[j] = indexReservoir[reservoirCount];
                continue;
            }
            j++;
        }
        treeIteratorPhase.setPhase(phase);
        treeIteratorPhase.reset();
        // loop over all the atoms.  Any atoms whose index is 
        for (IAtom a = treeIteratorPhase.nextAtom(); a != null;
             a = treeIteratorPhase.nextAtom()) {
            if (a.getGlobalIndex() > maxIndex-reservoirSize) {
                PhaseAtomIndexChangedEvent event = new PhaseAtomIndexChangedEvent(phase, a, a.getGlobalIndex());
                // Just re-invoke the Atom's method without first "returning"
                // the index to the reservoir.  The old index gets dropped on the
                // floor.
                a.setGlobalIndex(this);
                phaseEventManager.fireEvent(event);
            }
        }
        maxIndex -= reservoirSize;
        PhaseGlobalAtomIndexEvent event = new PhaseGlobalAtomIndexEvent(phase, maxIndex);
        phaseEventManager.fireEvent(event);
        if (reservoirCount != 0) {
            System.out.println("reservoir still has atoms:");
            for (int i=0; i<reservoirCount; i++) {
                System.out.print(indexReservoir[i]+" ");
            }
            throw new RuntimeException("I was fully expecting the reservoir to be empty!");
        }
    }

    /**
     * Notifies the SpeciesMaster that the given number of new Atoms will be
     * added to the system.  It's not required to call this method before
     * adding atoms, but if adding many Atoms, calling this will improve
     * performance.
     */
    public void notifyNewAtoms(int numNewAtoms) {
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // max index has actually increased, there's no harm since there's 
        // nothing that says the max index can't be too large.
        if (numNewAtoms > reservoirCount) {
            PhaseGlobalAtomIndexEvent event = new PhaseGlobalAtomIndexEvent(phase, maxIndex + numNewAtoms - reservoirCount);
            phaseEventManager.fireEvent(event);
        }            
    }
    
    /**
     * Sets the size of the atom global index reservoir.
     * @param size
     */
    public void setIndexReservoirSize(int size) {
        if (size < 0) {
            throw new IllegalArgumentException("Reservoir size must not be negative");
        }
        collapseGlobalIndices();
        // Set the actual reservoir size to one more because we collapse the
        // indices when it's full, not when it's full and we have another to add.
        reservoirSize = size+1;
        indexReservoir = new int[reservoirSize];
    }

    /**
     * Returns the size of the reservoir; the number of Atom that can be
     * removed without triggering an index collapse.
     */
    public int getIndexReservoirSize() {
        return reservoirSize-1;
    }

    public void addAtomNotify(IAtom newAtom) {
        if (newAtom.getParentGroup() instanceof SpeciesAgent) {
            moleculeCount++;
        } else if (newAtom instanceof SpeciesAgent) {
            moleculeCount += ((SpeciesAgent) newAtom)
                    .getNMolecules();
        }

        newAtom.setGlobalIndex(this);
        if (newAtom.isLeaf()) {
            ((AtomLeaf)newAtom).setLeafIndex(leafList.size());
            leafList.add(newAtom);
        } else {
            treeIteratorRoot.setRootAtom(newAtom);
            treeIteratorRoot.reset();
            for (IAtom childAtom = treeIteratorRoot.nextAtom(); childAtom != null;
                 childAtom = treeIteratorRoot.nextAtom()) {
                if (childAtom.getType().isLeaf()) {
                    ((AtomLeaf)childAtom).setLeafIndex(leafList.size());
                    leafList.add(childAtom);
                }
                childAtom.setGlobalIndex(this);
            }
        }
        phaseEventManager.fireEvent(new PhaseAtomAddedEvent(phase, newAtom));
    }

    //updating of leaf atomList may not be efficient enough for repeated
    // use, but is probably ok
    public void removeAtomNotify(IAtom oldAtom) {
        if (oldAtom.getParentGroup() instanceof SpeciesAgent) {
            moleculeCount--;
        } else if (oldAtom instanceof SpeciesAgent) {
            moleculeCount -= ((SpeciesAgent)oldAtom).getNMolecules();
//            ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
        }
        
        phaseEventManager.fireEvent(new PhaseAtomRemovedEvent(phase, oldAtom));
        returnGlobalIndex(oldAtom.getGlobalIndex());
        if (oldAtom.isLeaf()) {
            int leafIndex = ((AtomLeaf)oldAtom).getLeafIndex();
            leafList.removeAndReplace(leafIndex);
            leafList.maybeTrimToSize();
            // if we removed didn't remove the last atom, removeAndReplace
            // inserted the last atom in the emtpy spot.  Set its leaf index.
            if (leafList.size() > leafIndex) {
                ((AtomLeaf)leafList.get(leafIndex)).setLeafIndex(leafIndex);
            }
        } else {
            treeIteratorRoot.setRootAtom(oldAtom);
            treeIteratorRoot.reset();
            for (IAtom childAtom = treeIteratorRoot.nextAtom(); childAtom != null;
                 childAtom = treeIteratorRoot.nextAtom()) {
                returnGlobalIndex(childAtom.getGlobalIndex());
                if (childAtom.getType().isLeaf()) {
                    int leafIndex = ((AtomLeaf)childAtom).getLeafIndex();
                    leafList.removeAndReplace(leafIndex);
                    if (leafList.size() > leafIndex) {
                        ((AtomLeaf)leafList.get(leafIndex)).setLeafIndex(leafIndex);
                    }
                }
            }
            leafList.maybeTrimToSize();
        }
    }

    private static final long serialVersionUID = 2L;
    private final Phase phase;
    private final AtomArrayList agentList;
    /**
     * List of leaf atoms in phase
     */
    private final AtomArrayList leafList = new AtomArrayList();

    private final AtomIteratorTreePhase treeIteratorPhase;
    private final AtomIteratorTreeRoot treeIteratorRoot;

    protected int moleculeCount;
    protected final PhaseEventManager phaseEventManager;

    private int[] indexReservoir;
    private int reservoirSize = 50;
    private int reservoirCount;
    private int maxIndex;
    
    /**
     * non-graphic main method to test handling of leaf atom list.
     */
    public static void main(String args[]) {

        Simulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        sim.getSpeciesManager().addSpecies(species2);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species0);
        Phase phase = new Phase(sim);
        phase.getAgent(species0).setNMolecules(4);
        phase.getAgent(species1).setNMolecules(2);
        phase.getAgent(species2).setNMolecules(2);

        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
            System.out.println(a.toString());
        }
        System.out.println();
    }//end of main

}//end of SpeciesMaster
