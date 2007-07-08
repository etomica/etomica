package etomica.atom;

import etomica.atom.iterator.AtomIteratorTreeBox;
import etomica.atom.iterator.AtomIteratorTreeRoot;
import etomica.box.Box;
import etomica.box.BoxAtomAddedEvent;
import etomica.box.BoxAtomIndexChangedEvent;
import etomica.box.BoxAtomLeafIndexChangedEvent;
import etomica.box.BoxAtomRemovedEvent;
import etomica.box.BoxEventManager;
import etomica.box.BoxGlobalAtomIndexEvent;
import etomica.box.BoxGlobalAtomLeafIndexEvent;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * Coordinator of all species agents in a box. Parent is SpeciesRoot, and
 * children are SpeciesAgent instances. All instances of SpeciesMaster in a
 * given simulation share the same AtomType.
 * 
 * @author David Kofke and Andrew Schultz
 */

public final class AtomManager implements java.io.Serializable {

    public AtomManager(Box p, BoxEventManager eventManager) {
        agentList = new AtomArrayList();
        boxEventManager = eventManager;
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        leafIndices = new int[0];
        reservoirCount = 0;
        box = p;
        treeIteratorRoot = new AtomIteratorTreeRoot();
        treeIteratorRoot.setDoAllNodes(true);
        treeIteratorRoot.setIterationDepth(Integer.MAX_VALUE);
        treeIteratorBox = new AtomIteratorTreeBox();
    }

    /**
     * Adds the given SpeciesAgent to the SpeciesMaster's list of
     * SpeciesAgents.  This method should be called by Species.
     */
    public void addSpeciesAgent(ISpeciesAgent newSpeciesAgent) {
        newSpeciesAgent.setIndex(agentList.getAtomCount());
        agentList.add(newSpeciesAgent);
        addAtomNotify(newSpeciesAgent);
    }
    
    /**
     * Notifies the SpeciesMaster that a Species has been removed.  This method
     * should only be called by the SpeciesManager.
     */
    public void removeSpeciesNotify(Species species) {
        for (int i=0; i<agentList.getAtomCount(); i++) {
            IAtom speciesAgent = agentList.getAtom(i);
            if (speciesAgent.getType().getSpecies() == species) {
                agentList.removeAndReplace(i);
                if (agentList.getAtomCount() > i) {
                    agentList.getAtom(i).setIndex(i);
                }
                agentList.maybeTrimToSize();
                removeAtomNotify(speciesAgent);
            }
        }
    }

    public AtomSet getAgentList() {
        return agentList;
    }

    /**
     * Returns the number of molecules in the Box
     */
    public int moleculeCount() {
        return moleculeCount;
    }

    /**
     * Returns an AtomArrayList containing the leaf atoms in the Box
     */
    public AtomSet getLeafList() {
        return leafList;
    }
    
    public Box getBox() {
        return box;
    }

    /**
     * Returns a "global" index for the Box.  This method should only be
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
     * Returns the maximum global index for the Box.
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
        treeIteratorBox.setBox(box);
        treeIteratorBox.reset();
        // loop over all the atoms.  Any atoms whose index is larger than what
        // the new maxIndex will be get new indices
        for (IAtom a = treeIteratorBox.nextAtom(); a != null;
             a = treeIteratorBox.nextAtom()) {
            if (a.getGlobalIndex() > maxIndex-reservoirSize) {
                // Just re-invoke the Atom's method without first "returning"
                // the index to the reservoir.  The old index gets dropped on the
                // floor.
                int oldGlobalIndex = a.getGlobalIndex();
                BoxAtomIndexChangedEvent event = new BoxAtomIndexChangedEvent(box, a, oldGlobalIndex);
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("reassigning global index for "+a);
                }
                a.setGlobalIndex(this);
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("reassigned global index for "+a+" from "+oldGlobalIndex+" to "+a.getGlobalIndex());
                }
                leafIndices[a.getGlobalIndex()] = leafIndices[oldGlobalIndex];
                boxEventManager.fireEvent(event);
            }
        }
        maxIndex -= reservoirSize;
        if (leafIndices.length > maxIndex + 1 + reservoirSize) {
            leafIndices = Arrays.resizeArray(leafIndices, maxIndex+1);
        }
        BoxGlobalAtomIndexEvent event = new BoxGlobalAtomIndexEvent(box, maxIndex);
        boxEventManager.fireEvent(event);
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
    public void notifyNewAtoms(int numNewAtoms, int numNewLeafAtoms) {
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // actual max index has already increased, there's no harm since
        // there's nothing that says the max index can't be too large.
        if (numNewAtoms > reservoirCount) {
            BoxGlobalAtomIndexEvent event = new BoxGlobalAtomIndexEvent(box, maxIndex + numNewAtoms - reservoirCount);
            boxEventManager.fireEvent(event);
            leafIndices = Arrays.resizeArray(leafIndices, maxIndex + numNewAtoms - reservoirCount + 1 + reservoirSize);
        }
        if (numNewLeafAtoms > 1) {
            BoxGlobalAtomLeafIndexEvent leafEvent = new BoxGlobalAtomLeafIndexEvent(box, leafList.getAtomCount() + numNewLeafAtoms);
            boxEventManager.fireEvent(leafEvent);
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
        if (newAtom.getType().getDepth() == AtomAddressManager.MOLECULE_DEPTH) {
            moleculeCount++;
        } else if (newAtom instanceof ISpeciesAgent) {
            moleculeCount += ((ISpeciesAgent) newAtom)
                    .getNMolecules();
        }

        newAtom.setGlobalIndex(this);
        if (!(newAtom instanceof IAtomGroup)) {
            int globalIndex = newAtom.getGlobalIndex();
            if (globalIndex > leafIndices.length-1) {
                leafIndices = Arrays.resizeArray(leafIndices, globalIndex + 1 + reservoirSize);
            }
            leafIndices[globalIndex] = leafList.getAtomCount();
            leafList.add(newAtom);
        } else {
            treeIteratorRoot.setRootAtom(newAtom);
            treeIteratorRoot.reset();
            for (IAtom childAtom = treeIteratorRoot.nextAtom(); childAtom != null;
                 childAtom = treeIteratorRoot.nextAtom()) {
                childAtom.setGlobalIndex(this);
                if (childAtom.getType().isLeaf()) {
                    int globalIndex = childAtom.getGlobalIndex();
                    if (globalIndex > leafIndices.length-1) {
                        leafIndices = Arrays.resizeArray(leafIndices, globalIndex + 1 + reservoirSize);
                    }
                    leafIndices[globalIndex] = leafList.getAtomCount();
                    leafList.add(childAtom);
                }
            }
        }
        boxEventManager.fireEvent(new BoxAtomAddedEvent(box, newAtom));
    }

    //updating of leaf atomList may not be efficient enough for repeated
    // use, but is probably ok
    public void removeAtomNotify(IAtom oldAtom) {
        if (oldAtom.getType().getDepth() == AtomAddressManager.MOLECULE_DEPTH) {
            moleculeCount--;
        } else if (oldAtom instanceof ISpeciesAgent) {
            moleculeCount -= ((ISpeciesAgent)oldAtom).getNMolecules();
//            ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
        }
        
        boxEventManager.fireEvent(new BoxAtomRemovedEvent(box, oldAtom));
        if (!(oldAtom instanceof IAtomGroup)) {
            int leafIndex = leafIndices[oldAtom.getGlobalIndex()];
            leafList.removeAndReplace(leafIndex);
            leafList.maybeTrimToSize();
            // if we didn't remove the last atom, removeAndReplace
            // inserted the last atom in the emtpy spot.  Set its leaf index.
            if (leafList.getAtomCount() > leafIndex) {
                IAtom movedAtom = leafList.getAtom(leafIndex);
                int globalIndex = movedAtom.getGlobalIndex();
                BoxAtomLeafIndexChangedEvent event = new BoxAtomLeafIndexChangedEvent(box, movedAtom, leafIndices[globalIndex]);
                leafIndices[globalIndex] = leafIndex;
                boxEventManager.fireEvent(event);
            }
            returnGlobalIndex(oldAtom.getGlobalIndex());
        } else {
            returnGlobalIndex(oldAtom.getGlobalIndex());
            treeIteratorRoot.setRootAtom(oldAtom);
            treeIteratorRoot.reset();
            for (IAtom childAtom = treeIteratorRoot.nextAtom(); childAtom != null;
                 childAtom = treeIteratorRoot.nextAtom()) {
                if (childAtom.getType().isLeaf()) {
                    int leafIndex = leafIndices[childAtom.getGlobalIndex()];
                    leafList.removeAndReplace(leafIndex);
                    if (leafList.getAtomCount() > leafIndex) {
                        IAtom movedAtom = leafList.getAtom(leafIndex);
                        int globalIndex = movedAtom.getGlobalIndex();
                        BoxAtomLeafIndexChangedEvent event = new BoxAtomLeafIndexChangedEvent(box, movedAtom, leafIndices[globalIndex]);
                        leafIndices[globalIndex] = leafIndex;
                        boxEventManager.fireEvent(event);
                    }
                }
                returnGlobalIndex(childAtom.getGlobalIndex());
            }
            leafList.maybeTrimToSize();
        }
    }
    
    /**
     * Returns the index of the given leaf atom within the SpeciesMaster's
     * leaf list.  The given leaf atom must be in the SpeciesMaster's Box. 
     */
    public int getLeafIndex(IAtom atomLeaf) {
        return leafIndices[atomLeaf.getGlobalIndex()];
    }

    private static final long serialVersionUID = 2L;
    private final Box box;
    private final AtomArrayList agentList;
    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList = new AtomArrayList();

    protected final AtomIteratorTreeBox treeIteratorBox;
    protected final AtomIteratorTreeRoot treeIteratorRoot;

    protected int moleculeCount;
    protected final BoxEventManager boxEventManager;

    protected int[] indexReservoir;
    protected int reservoirSize = 50;
    protected int reservoirCount;
    protected int maxIndex;
    
    protected int[] leafIndices;
    
    /**
     * non-graphic main method to test handling of leaf atom list.
     */
    public static void main(String args[]) {

        ISimulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        sim.getSpeciesManager().addSpecies(species2);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species0);
        Box box = new Box(sim);
        sim.addBox(box);
        box.getAgent(species0).setNMolecules(4);
        box.getAgent(species1).setNMolecules(2);
        box.getAgent(species2).setNMolecules(2);

        AtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned a = (IAtomPositioned)leafList.getAtom(iLeaf);
            System.out.println(a.toString());
        }
        System.out.println();
    }

}
