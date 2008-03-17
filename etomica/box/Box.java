package etomica.box;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IBoxEventManager;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;

import etomica.action.BoxInflate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtomLeaf;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * A Box collects all atoms that interact with one another; atoms in different
 * boxs do not interact. These are the important features of a Box:
 * <p>
 * <ol>
 * <li>It holds a SpeciesMaster instance, which provides the root of a
 * hierarchy of atoms that represent the physical objects that interact.
 * <li>It holds a Boundary object, obtained from the governing Space, that
 * defines the volume of the box and the behavior of atoms as they move into or
 * across the boundary of the box.
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the box (such as a change in its boundary).
 * <li>Each Box has a unique index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Box
 * construction. This index is useful when collecting things in reference to the
 * box.
 * </ol>
 * A box is acted upon by an Integrator instance to move its atoms around and
 * generate configurations. Properties of a box are measured by MeterAbstract
 * instances which are simply DataSource objects that require a box to
 * generate their data. <br>
 * A simulation may involve more than one box. All Box instances are
 * registered with the simulation specified upon their construction, and
 * may be accessed via the simulation's getBoxs method.
 * 
 * @author David Kofke, Andrew Schultz
 * @see Boundary
 */
public class Box implements java.io.Serializable, IBox {

    /**
     * Constructs box with default rectangular periodic boundary.
     */
    public Box(ISimulation sim, Space _space) {
        this(new BoundaryRectangularPeriodic(sim.getRandom(), _space), _space);
    }
    
    /**
     * Constructs box with the given boundary
     */
    public Box(IBoundary boundary, Space _space) {
    	this.space = _space;
        eventManager = new BoxEventManager();
        setBoundary(boundary);
        
        inflateEvent = new BoxInflateEvent(this);
        moleculeLists = new AtomArrayList[0];
        allMoleculeList = new AtomSetAllMolecules();
        allMoleculeList.setMoleculeLists(moleculeLists);
        leafList = new AtomArrayList();
        
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        leafIndices = new int[0];
        reservoirCount = 0;
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#resetIndex(etomica.simulation.ISimulation)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#resetIndex(etomica.api.ISimulation)
	 */
    public void resetIndex(ISimulation sim) {
        if (sim == null) {
            // sim is notifying us that we got removed.
            index = 0;
            return;
        }
        IBox[] boxs = sim.getBoxs();
        for (int i=0; i<boxs.length; i++) {
            if (boxs[i] == this) {
                index = i;
                return;
            }
        }
        // you really shouldn't be calling resetIndex unless you're a simulation!
        throw new IllegalArgumentException(sim+" does not contain me");
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#getIndex()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getIndex()
	 */
    public int getIndex() {
        return index;
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#toString()
	 */
    public String toString() {
        return "Box"+getIndex();
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#addNewMolecule(etomica.species.ISpecies)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#addNewMolecule(etomica.api.ISpecies)
	 */
    public IMolecule addNewMolecule(ISpecies species) {
        IMolecule aNew = species.makeMolecule();
        addMolecule(aNew, species);
        return aNew;
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#addMolecule(etomica.atom.IMolecule)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#addMolecule(etomica.api.IMolecule)
	 */
    public void addMolecule(IMolecule molecule) {
        ISpecies species = molecule.getType().getSpecies();
        addMolecule(molecule, species);
    }
    
    protected void addMolecule(IMolecule molecule, ISpecies species) {
        int speciesIndex = species.getIndex();
        if (moleculeLists[speciesIndex].contains(molecule)) {
            throw new RuntimeException("you bastard!");
        }
        molecule.setIndex(moleculeLists[speciesIndex].getAtomCount());
        moleculeLists[speciesIndex].add(molecule);
        allMoleculeList.setMoleculeLists(moleculeLists);
        addAtomNotify(molecule);
        for (int i=0; i<moleculeLists[speciesIndex].getAtomCount(); i++) {
            if (moleculeLists[speciesIndex].getAtom(i).getIndex() != i) {
                throw new RuntimeException("oops "+molecule+" "+moleculeLists[speciesIndex].getAtom(i)+" "+i);
            }
        }
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#removeMolecule(etomica.atom.IMolecule)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#removeMolecule(etomica.api.IMolecule)
	 */
    public void removeMolecule(IMolecule molecule) {
        ISpecies species = molecule.getType().getSpecies();
        removeMolecule(molecule, species);
    }
    
    protected void removeMolecule(IMolecule molecule, ISpecies species) {
        int moleculeIndex = molecule.getIndex();
        AtomArrayList moleculeList = moleculeLists[species.getIndex()];
        if (moleculeList.getAtom(moleculeIndex) != molecule) {
            throw new IllegalArgumentException("can't find "+molecule);
        }
        if (moleculeIndex < moleculeList.getAtomCount()-1) {
            moleculeList.removeAndReplace(moleculeIndex);
            moleculeList.getAtom(moleculeIndex).setIndex(moleculeIndex);
        }
        else {
            moleculeList.remove(moleculeIndex);
        }
        allMoleculeList.setMoleculeLists(moleculeLists);
        removeAtomNotify(molecule);
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#setNMolecules(etomica.species.ISpecies, int)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#setNMolecules(etomica.api.ISpecies, int)
	 */
    public void setNMolecules(ISpecies species, int n) {
        int speciesIndex = species.getIndex();
        AtomArrayList moleculeList = moleculeLists[speciesIndex];
        int currentNMolecules = moleculeList.getAtomCount();
        notifyNewAtoms((n-currentNMolecules)*(1+species.getNumLeafAtoms()),
                                     (n-currentNMolecules)*species.getNumLeafAtoms());
        if(n < 0) {
            throw new IllegalArgumentException("Number of molecules cannot be negative");
        }
        for(int i=currentNMolecules; i<n; i++) {
            addNewMolecule(species);
        }
        for (int i=currentNMolecules; i>n; i--) {
            removeMolecule((IMolecule)moleculeList.getAtom(i-1));
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getNMolecules(etomica.species.ISpecies)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getNMolecules(etomica.api.ISpecies)
	 */
    public int getNMolecules(ISpecies species) {
        int speciesIndex = species.getIndex();
        return moleculeLists[speciesIndex].getAtomCount();
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMoleculeList(etomica.species.ISpecies)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMoleculeList(etomica.api.ISpecies)
	 */
    public IAtomSet getMoleculeList(ISpecies species) {
        return moleculeLists[species.getIndex()];
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMoleculeList()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMoleculeList()
	 */
    public IAtomSet getMoleculeList() {
        return allMoleculeList;
    }
  
    /* (non-Javadoc)
	 * @see etomica.box.IBox#setBoundary(etomica.space.Boundary)
	 */
     /* (non-Javadoc)
	 * @see etomica.box.IBox#setBoundary(etomica.api.IBoundary)
	 */
    public void setBoundary(IBoundary b) {
        boundary = b;
     }
     
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getBoundary()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getBoundary()
	 */
    public final IBoundary getBoundary() {return boundary;}
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#setDimensions(etomica.api.IVector)
	 */
    public final void setDimensions(IVector d) {
        boundary.setDimensions(d);
        eventManager.fireEvent(inflateEvent);
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#volume()
	 */
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#setDensity(double)
	 */
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        BoxInflate inflater = new BoxInflate(this, space);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getDensity()
	 */
    public double getDensity() {return moleculeCount()/boundary.volume();}
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getDensityDimension()
	 */
    public Dimension getDensityDimension() {
        return new DimensionRatio("Density",Quantity.DIMENSION,Volume.DIMENSION);
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#atomCount()
	 */
    public int atomCount() {return leafList.getAtomCount();}

    /* (non-Javadoc)
	 * @see etomica.box.IBox#getEventManager()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getEventManager()
	 */
    public IBoxEventManager getEventManager() {
        return eventManager;
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#addSpeciesNotify(etomica.api.ISpecies)
	 */
    public void addSpeciesNotify(ISpecies species) {
        moleculeLists = (AtomArrayList[])Arrays.addObject(moleculeLists, new AtomArrayList());
        allMoleculeList.setMoleculeLists(moleculeLists);
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#removeSpeciesNotify(etomica.api.ISpecies)
	 */
    public void removeSpeciesNotify(ISpecies species) {
        moleculeLists = (AtomArrayList[])Arrays.removeObject(moleculeLists, moleculeLists[species.getIndex()]);
        allMoleculeList.setMoleculeLists(moleculeLists);
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#moleculeCount()
	 */
    public int moleculeCount() {
        return allMoleculeList.getAtomCount();
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#getLeafList()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getLeafList()
	 */
    public IAtomSet getLeafList() {
        return leafList;
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#requestGlobalIndex()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#requestGlobalIndex()
	 */
    public int requestGlobalIndex() {
        if (reservoirCount == 0) {
            return ++maxIndex;
        }
        reservoirCount--;
        return indexReservoir[reservoirCount];
    }
    
    protected void returnGlobalIndex(int atomIndex) {
        // add the index to the reservoir first
        indexReservoir[reservoirCount] = atomIndex;
        reservoirCount++;
        if (reservoirCount == reservoirSize) {
            // reservoir is full
            collapseGlobalIndices();
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMaxGlobalIndex()
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getMaxGlobalIndex()
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
        if (Debug.ON) {
            int[] oldIndexArray = new int[maxIndex+1];
            for (int i=0; i<allMoleculeList.getAtomCount(); i++) {
                IAtom a = allMoleculeList.getAtom(i);
                if (oldIndexArray[a.getGlobalIndex()] != 0) {
                    throw new RuntimeException("molecule double index "+a.getGlobalIndex());
                }
                oldIndexArray[a.getGlobalIndex()] = 1;
            }
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom a = leafList.getAtom(i);
                if (oldIndexArray[a.getGlobalIndex()] != 0) {
                    throw new RuntimeException("leaf double index "+a.getGlobalIndex());
                }
                oldIndexArray[a.getGlobalIndex()] = 2;
            }
            for (int i=0; i<reservoirCount; i++) {
                if (oldIndexArray[indexReservoir[i]] != 0) {
                    throw new RuntimeException("reservoir double index "+indexReservoir[i]);
                }
                oldIndexArray[indexReservoir[i]] = 3;
            }
            for (int i=0; i<maxIndex+1; i++) {
                if (oldIndexArray[i] == 0) {
                    throw new RuntimeException("where is "+i);
                }
            }
        }
        for (int j=0; j<reservoirCount; ) {
            if (indexReservoir[j] > maxIndex-reservoirSize) {
                // this index isn't useful to us, so just drop it
                reservoirCount--;
                indexReservoir[j] = indexReservoir[reservoirCount];
                continue;
            }
            j++;
        }

        // loop over all the atoms.  Any atoms whose index is larger than what
        // the new maxIndex will be get new indices
        int recycledCount = 0;
        // we can't loop over the actual atoms in the box because we might be
        // in the middle of removing a molecule.  the leaf list and molecule list
        // will contain the atoms that haven't had their global indices recycled.
        boolean isLeafList = allMoleculeList.getAtomCount() == 0;
        IAtomSet currentList = isLeafList ? leafList : allMoleculeList;
        for (int j = 0; j<currentList.getAtomCount(); j++) {
            IAtom a = currentList.getAtom(j);
            if (a.getGlobalIndex() > maxIndex-reservoirSize) {
                recycledCount++;
                // Just re-invoke the Atom's method without first "returning"
                // the index to the reservoir.  The old index gets dropped on the
                // floor.
                int oldGlobalIndex = a.getGlobalIndex();
                BoxAtomIndexChangedEvent event = new BoxAtomIndexChangedEvent(this, a, oldGlobalIndex);
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("reassigning global index for "+a);
                }
                a.setGlobalIndex(this);
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("reassigned global index for "+a+" from "+oldGlobalIndex+" to "+a.getGlobalIndex());
                }
                leafIndices[a.getGlobalIndex()] = leafIndices[oldGlobalIndex];
                eventManager.fireEvent(event);
            }
            else if (Debug.ON) {
                for (int i=0; i<reservoirCount; i++) {
                    if (a.getGlobalIndex() == indexReservoir[i]) {
                        throw new RuntimeException(a+" "+a.getGlobalIndex()+" was in the reservoir");
                    }
                }
            }
            if (!isLeafList && j == currentList.getAtomCount()-1) {
                // wow, this is evil.  we finished with the molecules, so now
                // go through the leaf atoms
                isLeafList = true;
                currentList = leafList;
                j = -1;
            }
        }
        maxIndex -= reservoirSize;
        if (leafIndices.length > maxIndex + 1 + reservoirSize) {
            leafIndices = Arrays.resizeArray(leafIndices, maxIndex+1);
        }
        BoxGlobalAtomIndexEvent event = new BoxGlobalAtomIndexEvent(this, maxIndex);
        eventManager.fireEvent(event);
        if (reservoirCount != 0) {
            System.out.println("reservoir still has atoms:");
            for (int i=0; i<reservoirCount; i++) {
                System.out.print(indexReservoir[i]+" ");
            }
            System.out.println();
            for (int i=0; i<leafList.getAtomCount(); i++) {
                System.out.println(leafList.getAtom(i)+" "+leafList.getAtom(i).getGlobalIndex());
            }
            for (int i=0; i<allMoleculeList.getAtomCount(); i++) {
                System.out.println(allMoleculeList.getAtom(i)+" "+allMoleculeList.getAtom(i).getGlobalIndex());
            }
            for (int i=0; i<reservoirCount; i++) {
                System.out.println("reservoir_"+i+" "+indexReservoir[i]);
            }
            throw new RuntimeException("I was fully expecting the reservoir to be empty!");
        }
    }

    /* (non-Javadoc)
	 * @see etomica.box.IBox#notifyNewAtoms(int, int)
	 */
    public void notifyNewAtoms(int numNewAtoms, int numNewLeafAtoms) {
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // actual max index has already increased, there's no harm since
        // there's nothing that says the max index can't be too large.
        if (numNewAtoms > reservoirCount) {
            BoxGlobalAtomIndexEvent event = new BoxGlobalAtomIndexEvent(this, maxIndex + numNewAtoms - reservoirCount);
            eventManager.fireEvent(event);
            leafIndices = Arrays.resizeArray(leafIndices, maxIndex + numNewAtoms - reservoirCount + 1 + reservoirSize);
        }
        if (numNewLeafAtoms > 1) {
            BoxGlobalAtomLeafIndexEvent leafEvent = new BoxGlobalAtomLeafIndexEvent(this, leafList.getAtomCount() + numNewLeafAtoms);
            eventManager.fireEvent(leafEvent);
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

    /* (non-Javadoc)
	 * @see etomica.box.IBox#addAtomNotify(etomica.api.IAtom)
	 */
    public void addAtomNotify(IAtom newAtom) {
        newAtom.setGlobalIndex(this);
        if (newAtom instanceof IAtomLeaf) {
            int globalIndex = newAtom.getGlobalIndex();
            if (globalIndex > leafIndices.length-1) {
                leafIndices = Arrays.resizeArray(leafIndices, globalIndex + 1 + reservoirSize);
            }
            leafIndices[globalIndex] = leafList.getAtomCount();
            leafList.add(newAtom);
        } else {
            IAtomSet childList = ((IMolecule)newAtom).getChildList();
            for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
                IAtomLeaf childAtom = (IAtomLeaf)childList.getAtom(iChild);
                childAtom.setGlobalIndex(this);
                int globalIndex = childAtom.getGlobalIndex();
                if (globalIndex > leafIndices.length-1) {
                    leafIndices = Arrays.resizeArray(leafIndices, globalIndex + 1 + reservoirSize);
                }
                leafIndices[globalIndex] = leafList.getAtomCount();
                leafList.add(childAtom);
            }
        }
        eventManager.fireEvent(new BoxAtomAddedEvent(this, newAtom));
    }

    //updating of leaf atomList may not be efficient enough for repeated
    // use, but is probably ok
    /* (non-Javadoc)
	 * @see etomica.box.IBox#removeAtomNotify(etomica.api.IAtom)
	 */
    public void removeAtomNotify(IAtom oldAtom) {
        eventManager.fireEvent(new BoxAtomRemovedEvent(this, oldAtom));
        if (oldAtom instanceof IAtomLeaf) {
            int leafIndex = leafIndices[oldAtom.getGlobalIndex()];
            leafList.removeAndReplace(leafIndex);
            leafList.maybeTrimToSize();
            // if we didn't remove the last atom, removeAndReplace
            // inserted the last atom in the emtpy spot.  Set its leaf index.
            if (leafList.getAtomCount() > leafIndex) {
                IAtom movedAtom = leafList.getAtom(leafIndex);
                int globalIndex = movedAtom.getGlobalIndex();
                BoxAtomLeafIndexChangedEvent event = new BoxAtomLeafIndexChangedEvent(this, movedAtom, leafIndices[globalIndex]);
                leafIndices[globalIndex] = leafIndex;
                eventManager.fireEvent(event);
            }
            returnGlobalIndex(oldAtom.getGlobalIndex());
        } else {
            IAtomSet childList = ((IMolecule)oldAtom).getChildList();
            returnGlobalIndex(oldAtom.getGlobalIndex());
            for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
                IAtomLeaf childAtom = (IAtomLeaf)childList.getAtom(iChild);
                int leafIndex = leafIndices[childAtom.getGlobalIndex()];
                leafList.removeAndReplace(leafIndex);
                if (leafList.getAtomCount() > leafIndex) {
                    IAtom movedAtom = leafList.getAtom(leafIndex);
                    int globalIndex = movedAtom.getGlobalIndex();
                    BoxAtomLeafIndexChangedEvent event = new BoxAtomLeafIndexChangedEvent(this, movedAtom, leafIndices[globalIndex]);
                    leafIndices[globalIndex] = leafIndex;
                    eventManager.fireEvent(event);
                }
                returnGlobalIndex(childAtom.getGlobalIndex());
            }
            leafList.maybeTrimToSize();
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getLeafIndex(etomica.atom.IAtom)
	 */
    /* (non-Javadoc)
	 * @see etomica.box.IBox#getLeafIndex(etomica.api.IAtom)
	 */
    public int getLeafIndex(IAtom atomLeaf) {
        return leafIndices[atomLeaf.getGlobalIndex()];
    }

    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList;

    protected int[] indexReservoir;
    protected int reservoirSize = 50;
    protected int reservoirCount;
    protected int maxIndex;
    
    protected int[] leafIndices;
    
    private static final long serialVersionUID = 2L;
    private IBoundary boundary;;
    private final IBoxEventManager eventManager;
    private final BoxEvent inflateEvent;
    protected AtomArrayList[] moleculeLists;
    private int index;
    private final Space space;
    protected final AtomSetAllMolecules allMoleculeList;
}
