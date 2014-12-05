/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.action.BoxInflate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IBoxEventManager;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculeArrayList;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
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
    public Box(ISpace _space) {
        this(new BoundaryRectangularPeriodic(_space), _space);
    }
    
    /**
     * Constructs box with the given boundary
     */
    public Box(IBoundary boundary, ISpace _space) {
    	this.space = _space;
        eventManager = new BoxEventManager(this);
        setBoundary(boundary);
        
        moleculeLists = new MoleculeArrayList[0];
        allMoleculeList = new AtomSetAllMolecules();
        allMoleculeList.setMoleculeLists(moleculeLists);
        leafList = new AtomArrayList();
        
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        reservoirCount = 0;
    }
    
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    public int getIndex() {
        return index;
    }
    
    public String toString() {
        return "Box"+getIndex();
    }
    
    public IMolecule addNewMolecule(ISpecies species) {
        IMolecule aNew = species.makeMolecule();
        addMolecule(aNew);
        return aNew;
    }

    public void addMolecule(IMolecule molecule) {
        int speciesIndex = molecule.getType().getIndex();
        if (Debug.ON) {
            for (int i=0; i<moleculeLists[speciesIndex].getMoleculeCount(); i++) {
                if (moleculeLists[speciesIndex].getMolecule(i) == molecule) {
                    throw new RuntimeException("you bastard!");
                }
            }
        }
        molecule.setIndex(moleculeLists[speciesIndex].getMoleculeCount());
        moleculeLists[speciesIndex].add(molecule);
        allMoleculeList.setMoleculeLists(moleculeLists);

        IAtomList childList = molecule.getChildList();
        int nLeafAtoms = leafList.getAtomCount();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtom childAtom = childList.getAtom(iChild);
            childAtom.setLeafIndex(nLeafAtoms++);
            leafList.add(childAtom);
        }
        eventManager.moleculeAdded(molecule);

        if (Debug.ON) {
            for (int i=0; i<moleculeLists[speciesIndex].getMoleculeCount(); i++) {
                if (moleculeLists[speciesIndex].getMolecule(i).getIndex() != i) {
                    throw new RuntimeException("oops "+molecule+" "+moleculeLists[speciesIndex].getMolecule(i)+" "+i);
                }
            }
        }
    }

    public void removeMolecule(IMolecule molecule) {
        int moleculeIndex = molecule.getIndex();
        MoleculeArrayList moleculeList = moleculeLists[molecule.getType().getIndex()];
        if (Debug.ON && moleculeList.getMolecule(moleculeIndex) != molecule) {
            throw new IllegalArgumentException("can't find "+molecule);
        }
        if (moleculeIndex < moleculeList.getMoleculeCount()-1) {
            moleculeList.removeAndReplace(moleculeIndex);
            IMolecule replacingMolecule = moleculeList.getMolecule(moleculeIndex);
            replacingMolecule.setIndex(moleculeIndex);
            eventManager.moleculeIndexChanged(replacingMolecule, moleculeList.getMoleculeCount());
        }
        else {
            moleculeList.remove(moleculeIndex);
        }
        allMoleculeList.setMoleculeLists(moleculeLists);

        eventManager.moleculeRemoved(molecule);
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtom childAtom = childList.getAtom(iChild);
            int leafIndex = childAtom.getLeafIndex();
            leafList.removeAndReplace(leafIndex);
            if (leafList.getAtomCount() > leafIndex) {
                IAtom movedAtom = leafList.getAtom(leafIndex);
                int movedLeafIndex = movedAtom.getLeafIndex();
                movedAtom.setLeafIndex(leafIndex);
                eventManager.atomLeafIndexChanged(movedAtom, movedLeafIndex);
            }
        }
        leafList.maybeTrimToSize();
    }
    
    public void setNMolecules(ISpecies species, int n) {
        int speciesIndex = species.getIndex();
        MoleculeArrayList moleculeList = moleculeLists[speciesIndex];
        int currentNMolecules = moleculeList.getMoleculeCount();
        int moleculeLeafAtoms = 0;
        IMolecule newMolecule0 = null;
        if (currentNMolecules > 0) {
            moleculeLeafAtoms = moleculeList.getMolecule(0).getChildList().getAtomCount();
        }
        else if (n > currentNMolecules) {
            newMolecule0 = species.makeMolecule();
            moleculeLeafAtoms = newMolecule0.getChildList().getAtomCount();
        }
        notifyNewMolecules(species, (n-currentNMolecules), moleculeLeafAtoms);
        if(n < 0) {
            throw new IllegalArgumentException("Number of molecules cannot be negative");
        }
        if (n > currentNMolecules) {
            moleculeLists[species.getIndex()].ensureCapacity(n);
            leafList.ensureCapacity(leafList.getAtomCount()+(n-currentNMolecules)*moleculeLeafAtoms);
            if (newMolecule0 != null) {
                addMolecule(newMolecule0);
                currentNMolecules++;
            }
            for(int i=currentNMolecules; i<n; i++) {
                addMolecule(species.makeMolecule());
            }
        }
        else {
            for (int i=currentNMolecules; i>n; i--) {
                removeMolecule(moleculeList.getMolecule(i-1));
            }
        }
    }
    
    public int getNMolecules(ISpecies species) {
        int speciesIndex = species.getIndex();
        return moleculeLists[speciesIndex].getMoleculeCount();
    }
    
    public IMoleculeList getMoleculeList(ISpecies species) {
        return moleculeLists[species.getIndex()];
    }

    public IMoleculeList getMoleculeList() {
        return allMoleculeList;
    }
  
    public void setBoundary(IBoundary b) {
        boundary = b;
        boundary.setBox(this);
     }
     
    public final IBoundary getBoundary() {return boundary;}
    
    public void setDensity(double rho) {
        double vNew = getMoleculeList().getMoleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        BoxInflate inflater = new BoxInflate(this, space);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }

    public IBoxEventManager getEventManager() {
        return eventManager;
    }

    public void addSpeciesNotify(ISpecies species) {
        moleculeLists = (MoleculeArrayList[])Arrays.addObject(moleculeLists, new MoleculeArrayList());
        allMoleculeList.setMoleculeLists(moleculeLists);
    }
    
    public void removeSpeciesNotify(ISpecies species) {
        moleculeLists = (MoleculeArrayList[])Arrays.removeObject(moleculeLists, moleculeLists[species.getIndex()]);
        allMoleculeList.setMoleculeLists(moleculeLists);
    }

    public IAtomList getLeafList() {
        return leafList;
    }
    
    /**
     * Notifies the SpeciesMaster that the given number of new Atoms will be
     * added to the system.  It's not required to call this method before
     * adding atoms, but if adding many Atoms, calling this will improve
     * performance.
     */
    protected void notifyNewMolecules(ISpecies species, int numNewMolecules, int moleculeLeafAtoms) {
        if (numNewMolecules < 1) return;
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // actual max index has already increased, there's no harm since
        // there's nothing that says the max index can't be too large.
        int numNewLeafAtoms = numNewMolecules * moleculeLeafAtoms;
        if (numNewLeafAtoms + numNewMolecules > reservoirCount) {
            eventManager.globalAtomIndexChanged(maxIndex + numNewMolecules + numNewLeafAtoms - reservoirCount);
        }
        eventManager.numberMolecules(species, moleculeLists[species.getIndex()].getMoleculeCount() + numNewMolecules);
        if (numNewLeafAtoms > 1) {
            eventManager.globalAtomLeafIndexChanged(leafList.getAtomCount() + numNewLeafAtoms);
        }
    }

    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList;

    protected int[] indexReservoir;
    protected int reservoirSize = 50;
    protected int reservoirCount;
    protected int maxIndex;
    
    private static final long serialVersionUID = 2L;
    private IBoundary boundary;;
    private final BoxEventManager eventManager;
    protected MoleculeArrayList[] moleculeLists;
    private int index;
    private final ISpace space;
    protected final AtomSetAllMolecules allMoleculeList;
}
