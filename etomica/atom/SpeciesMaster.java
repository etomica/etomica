package etomica.atom;

import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.phase.Phase;
import etomica.phase.PhaseEvent;
import etomica.phase.PhaseEventManager;
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

public final class SpeciesMaster extends Atom {

    protected int moleculeCount;
    public final static int SPECIES_TAB = AtomLinker.Tab.requestTabType();
    private final PhaseEventManager phaseEventManager;

    private final PhaseEvent changeIndexEvent;
    private final PhaseEvent maxGlobalIndexEvent;
    private int[] indexReservoir;
    private int reservoirSize = 50;
    private int reservoirCount;
    private int maxIndex;

    //reference to phase is kept in node

    /**
     * Tabbed list of leaf atoms in phase, suitable for iteration via an
     * AtomIteratorTabbedList.
     */
    public final AtomArrayList leafList = new AtomArrayList();

    public SpeciesMaster(Simulation sim, Phase p, PhaseEventManager eventManager) {
        super(sim.speciesRoot.getChildType(), new NodeFactory(p));
        phaseEventManager = eventManager;
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        reservoirCount = 0;
        setGlobalIndex(this);
        changeIndexEvent = new PhaseEvent(this,PhaseEvent.ATOM_CHANGE_INDEX);
        maxGlobalIndexEvent = new PhaseEvent(this,PhaseEvent.GLOBAL_INDEX);
    }

    //    public int atomCount() {return atomList.size();}//or could use
    // node.leafAtomCount()
    public int moleculeCount() {
        return moleculeCount;
    }

    public String signature() {
        Phase phase = node.parentPhase();
        return (phase != null) ? phase.getName()
                : "SpeciesMaster without phase";
    }
    
    public void removeSpecies(Species species) {
        AtomArrayList speciesAgents = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesAgents);
        iterator.reset();
        while (iterator.hasNext()) {
            Atom speciesAgent = iterator.nextAtom();
            if (speciesAgent.type.getSpecies() == species) {
                speciesAgent.node.dispose();
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
        AtomIteratorTree treeIterator = new AtomIteratorTree(this,Integer.MAX_VALUE,true);
        treeIterator.reset();
        // loop over all the atoms.  Any atoms whose index is 
        while (treeIterator.hasNext()) {
            Atom a = treeIterator.nextAtom();
            if (a.getGlobalIndex() > maxIndex-reservoirSize) {
                changeIndexEvent.setPhase(a.node.parentPhase());
                changeIndexEvent.setAtom(a);
                changeIndexEvent.setIndex(a.getGlobalIndex());
                // Just re-invoke the Atom's method without first "returning"
                // the index to the reservoir.  The old index gets dropped on the
                // floor.
                a.setGlobalIndex(this);
                phaseEventManager.fireEvent(changeIndexEvent);
            }
        }
        maxIndex -= reservoirSize;
        maxGlobalIndexEvent.setIndex(maxIndex);
        phaseEventManager.fireEvent(maxGlobalIndexEvent);
        if (reservoirCount != 0) {
            System.out.println("reservoir still has atoms:");
            for (int i=0; i<reservoirCount; i++) {
                System.out.print(indexReservoir[i]+" ");
            }
            throw new RuntimeException("I was fully expecting the reservoir to be empty!");
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

    private static final class MasterAtomTreeNode extends AtomTreeNodeGroup {

        MasterAtomTreeNode(Phase parentPhase, Atom atom) {
            super(atom);
            speciesMaster = (SpeciesMaster) atom;
            this.parentPhase = parentPhase;
            treeIterator.setDoAllNodes(true);
            treeIterator.setIterationDepth(Integer.MAX_VALUE);
            additionEvent = new PhaseEvent(atom, PhaseEvent.ATOM_ADDED);
            additionEvent.setPhase(parentPhase);
            removalEvent = new PhaseEvent(atom, PhaseEvent.ATOM_REMOVED);
            removalEvent.setPhase(parentPhase);
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

        /**
         * Returns true, because children are SpeciesAgent instances.
         */
        public final boolean childrenAreGroups() {
            return true;
        }

        public void addAtomNotify(Atom newAtom) {
            if (newAtom.node.parentGroup() instanceof SpeciesAgent) {
                speciesMaster.moleculeCount++;
            } else if (newAtom instanceof SpeciesAgent) {
                speciesMaster.moleculeCount += ((SpeciesAgent) newAtom)
                        .getNMolecules();
            }

            leafAtomCount += newAtom.node.leafAtomCount();
            if (newAtom.node.isLeaf()) {
                newAtom.setGlobalIndex((SpeciesMaster)atom);
                ((AtomTreeNodeLeaf)newAtom.node).setLeafIndex(speciesMaster.leafList.size());
                speciesMaster.leafList.add(newAtom);
            } else {
                treeIterator.setRoot(newAtom);
                treeIterator.reset();
                while (treeIterator.hasNext()) {
                    Atom childAtom = treeIterator.nextAtom();
                    if (childAtom.type.isLeaf()) {
                        ((AtomTreeNodeLeaf)childAtom.node).setLeafIndex(speciesMaster.leafList.size());
                        speciesMaster.leafList.add(childAtom);
                    }
                    childAtom.setGlobalIndex((SpeciesMaster)atom);
                }
            }
            additionEvent.setAtom(newAtom);
            ((SpeciesMaster)atom).phaseEventManager.fireEvent(additionEvent);
            if (parentNode() != null) {
                parentNode().addAtomNotify(newAtom);
            }
       }

        //updating of leaf atomList may not be efficient enough for repeated
        // use, but is probably ok
        public void removeAtomNotify(Atom oldAtom) {
            if (oldAtom.node.parentGroup() instanceof SpeciesAgent) {
                speciesMaster.moleculeCount--;
            } else if (oldAtom instanceof SpeciesAgent) {
                speciesMaster.moleculeCount -= ((SpeciesAgent) oldAtom)
                        .getNMolecules();
//                ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
            }
            
            if (oldAtom.node.isLeaf()) {
                leafAtomCount--;
                int leafIndex = ((AtomTreeNodeLeaf)oldAtom.node).getLeafIndex();
                speciesMaster.leafList.removeAndReplace(leafIndex);
                if (speciesMaster.leafList.size() > leafIndex) {
                    ((AtomTreeNodeLeaf)speciesMaster.leafList.get(leafIndex).node).setLeafIndex(leafIndex);
                }
            } else {
                leafAtomCount -= oldAtom.node.leafAtomCount();
                treeIterator.setRoot(oldAtom);
                treeIterator.reset();
                while (treeIterator.hasNext()) {
                    Atom childAtom = treeIterator.nextAtom();
                    if (childAtom.type.isLeaf()) {
                        int leafIndex = ((AtomTreeNodeLeaf)childAtom.node).getLeafIndex();
                        speciesMaster.leafList.removeAndReplace(leafIndex);
                        if (speciesMaster.leafList.size() > leafIndex) {
                            ((AtomTreeNodeLeaf)speciesMaster.leafList.get(leafIndex).node).setLeafIndex(leafIndex);
                        }
                    }
                }
            }
            removalEvent.setAtom(oldAtom);
            ((SpeciesMaster)atom).phaseEventManager.fireEvent(removalEvent);
            ((SpeciesMaster)atom).returnGlobalIndex(oldAtom.getGlobalIndex());
            if (parentNode() != null) {
                parentNode().removeAtomNotify(oldAtom);
            }
        }

        private final PhaseEvent additionEvent;
        private final PhaseEvent removalEvent;
        private final Phase parentPhase;
        private final SpeciesMaster speciesMaster;
        private final AtomIteratorTree treeIterator = new AtomIteratorTree();
    } //end of MasterAtomTreeNode

    private static final class NodeFactory implements AtomTreeNodeFactory, java.io.Serializable {

        Phase phase;

        NodeFactory(Phase p) {
            phase = p;
        }

        public AtomTreeNode makeNode(Atom atom) {
            return new MasterAtomTreeNode(phase, atom);
        }
    }

    /**
     * non-graphic main method to test handling of leaf atom list.
     */
    public static void main(String args[]) {

        Simulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        species0.setNMolecules(4);
        species1.setNMolecules(2);
        species2.setNMolecules(2);
        Phase phase = new Phase(sim);
        //        sim.elementCoordinator.go();

        AtomIteratorLeafAtoms leafIterator = new AtomIteratorLeafAtoms();
        leafIterator.setPhase(phase);
        leafIterator.reset();
        while (leafIterator.hasNext())
            System.out.println(leafIterator.next().toString());
        System.out.println();
    }//end of main

}//end of SpeciesMaster
