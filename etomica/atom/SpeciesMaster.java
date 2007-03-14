package etomica.atom;

import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorTree;
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

public final class SpeciesMaster extends AtomGroup {

    //reference to phase is kept in node

    /**
     * Tabbed list of leaf atoms in phase, suitable for iteration via an
     * AtomIteratorTabbedList.
     */
    public final AtomArrayList leafList = new AtomArrayList();

    public SpeciesMaster(Simulation sim, Phase p, PhaseEventManager eventManager) {
        super(sim.getSpeciesRoot().getSpeciesMasterType());
        phaseEventManager = eventManager;
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        reservoirCount = 0;
        setGlobalIndex(this);
        parentPhase = p;
        treeIterator.setDoAllNodes(true);
        treeIterator.setIterationDepth(Integer.MAX_VALUE);
    }

    //    public int atomCount() {return atomList.size();}//or could use
    // node.leafAtomCount()
    public int moleculeCount() {
        return moleculeCount;
    }

    public String signature() {
        Phase phase = parentPhase();
        return (phase != null) ? phase.getName()
                : "SpeciesMaster without phase";
    }
    
    public void removeSpecies(Species species) {
        AtomArrayList speciesAgents = getChildList();
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesAgents);
        iterator.reset();
        while (iterator.hasNext()) {
            Atom speciesAgent = iterator.nextAtom();
            if (speciesAgent.getType().getSpecies() == species) {
                speciesAgent.dispose();
            }
        }
    }

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
        treeIterator.setDoAllNodes(true);
        treeIterator.setRootAtom(this);
        treeIterator.reset();
        // loop over all the atoms.  Any atoms whose index is 
        while (treeIterator.hasNext()) {
            Atom a = treeIterator.nextAtom();
            if (a.getGlobalIndex() > maxIndex-reservoirSize) {
                PhaseAtomIndexChangedEvent event = new PhaseAtomIndexChangedEvent(parentPhase(), a, a.getGlobalIndex());
                // Just re-invoke the Atom's method without first "returning"
                // the index to the reservoir.  The old index gets dropped on the
                // floor.
                a.setGlobalIndex(this);
                phaseEventManager.fireEvent(event);
            }
        }
        maxIndex -= reservoirSize;
        PhaseGlobalAtomIndexEvent event = new PhaseGlobalAtomIndexEvent(parentPhase(), maxIndex);
        phaseEventManager.fireEvent(event);
        if (reservoirCount != 0) {
            System.out.println("reservoir still has atoms:");
            for (int i=0; i<reservoirCount; i++) {
                System.out.print(indexReservoir[i]+" ");
            }
            throw new RuntimeException("I was fully expecting the reservoir to be empty!");
        }
    }
    
    public void notifyNewAtoms(int numNewAtoms) {
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // max index has actually increased, there's no harm since there's 
        // nothing that says the max index can't be too large.
        if (numNewAtoms > reservoirCount) {
            PhaseGlobalAtomIndexEvent event = new PhaseGlobalAtomIndexEvent(parentPhase(), maxIndex + numNewAtoms - reservoirCount);
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
    
    public int getIndexReservoirSize() {
        return reservoirSize-1;
    }

    public Phase parentPhase() {
        return parentPhase;
    }

    public Species parentSpecies() {
        throw new RuntimeException(
                "Error:  Unexpected call to parentSpecies in SpeciesMaster");
    }

    public SpeciesAgent parentSpeciesAgent() {
        throw new RuntimeException(
                "Error:  Unexpected call to parentSpeciesAgent in SpeciesMaster");
    }

    /**
     * Throws a RuntimeException, because a species master is not contained
     * within a molecule.
     */
    public final Atom parentMolecule() {
        throw new RuntimeException(
                "Error:  Unexpected call to parentMolecule in SpeciesMaster");
    }

    public void addAtomNotify(Atom newAtom) {
        if (newAtom.parentGroup() instanceof SpeciesAgent) {
            moleculeCount++;
        } else if (newAtom instanceof SpeciesAgent) {
            moleculeCount += ((SpeciesAgent) newAtom)
                    .getNMolecules();
        }

        leafAtomCount += newAtom.leafAtomCount();
        if (newAtom.isLeaf()) {
            newAtom.setGlobalIndex(this);
            ((AtomLeaf)newAtom).setLeafIndex(leafList.size());
            leafList.add(newAtom);
        } else {
            treeIterator.setRootAtom(newAtom);
            treeIterator.reset();
            while (treeIterator.hasNext()) {
                Atom childAtom = treeIterator.nextAtom();
                if (childAtom.getType().isLeaf()) {
                    ((AtomLeaf)childAtom).setLeafIndex(leafList.size());
                    leafList.add(childAtom);
                }
                childAtom.setGlobalIndex(this);
            }
        }
        phaseEventManager.fireEvent(new PhaseAtomAddedEvent(parentPhase, newAtom));
        if (parent != null) {
            parent.addAtomNotify(newAtom);
        }
   }

    //updating of leaf atomList may not be efficient enough for repeated
    // use, but is probably ok
    public void removeAtomNotify(Atom oldAtom) {
        if (oldAtom.parentGroup() instanceof SpeciesAgent) {
            moleculeCount--;
        } else if (oldAtom instanceof SpeciesAgent) {
            moleculeCount -= ((SpeciesAgent)oldAtom).getNMolecules();
//            ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
        }
        
        phaseEventManager.fireEvent(new PhaseAtomRemovedEvent(parentPhase, oldAtom));
        if (oldAtom.isLeaf()) {
            leafAtomCount--;
            int leafIndex = ((AtomLeaf)oldAtom).getLeafIndex();
            returnGlobalIndex(oldAtom.getGlobalIndex());
            leafList.removeAndReplace(leafIndex);
            leafList.maybeTrimToSize();
            if (leafList.size() > leafIndex) {
                ((AtomLeaf)leafList.get(leafIndex)).setLeafIndex(leafIndex);
            }
        } else {
            leafAtomCount -= oldAtom.leafAtomCount();
            treeIterator.setRootAtom(oldAtom);
            treeIterator.reset();
            while (treeIterator.hasNext()) {
                Atom childAtom = treeIterator.nextAtom();
                returnGlobalIndex(childAtom.getGlobalIndex());
                if (childAtom.getType().isLeaf()) {
                    int leafIndex = ((AtomLeaf)childAtom).getLeafIndex();
                    leafList.removeAndReplace(leafIndex);
                    if (leafList.size() > leafIndex) {
                        ((AtomLeaf)leafList.get(leafIndex)).setLeafIndex(leafIndex);
                    }
                }
            }
        }
        if (parent != null) {
            parent.removeAtomNotify(oldAtom);
        }
    }


    private static final long serialVersionUID = 2L;
    private final Phase parentPhase;
    private final AtomIteratorTree treeIterator = new AtomIteratorTree();

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
        sim.getSpeciesRoot().addSpecies(species2);
        sim.getSpeciesRoot().addSpecies(species1);
        sim.getSpeciesRoot().addSpecies(species0);
        Phase phase = new Phase(sim);
        phase.getAgent(species0).setNMolecules(4);
        phase.getAgent(species1).setNMolecules(2);
        phase.getAgent(species2).setNMolecules(2);

        AtomIteratorLeafAtoms leafIterator = new AtomIteratorLeafAtoms();
        leafIterator.setPhase(phase);
        leafIterator.reset();
        while (leafIterator.hasNext())
            System.out.println(leafIterator.next().toString());
        System.out.println();
    }//end of main

}//end of SpeciesMaster
