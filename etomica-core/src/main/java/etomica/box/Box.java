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
 * <li>It holds lists of the molecules in the box, including a list of all
 * molecules, and an array of lists that can be selected by Species.</li>
 * <li>It holds a list of all atoms in the box.</li>
 * <li>It holds a Boundary object, obtained from the governing Space, that
 * defines the volume of the box and the behavior of atoms as they move into or
 * across the boundary of the box.</li>
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the box (such as a change in its boundary).</li>
 * <li>Each Box has a unique index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Box
 * construction. This index is useful when collecting things in reference to the
 * box.</li>
 * </ol>
 * The box maintains index values for each molecule and each atom. These values
 * are held as integer fields by the corresponding Molecule and Atom instances.
 * The indices are assigned when molecule (and its atoms) are added to the box via
 * the addMolecule method. The indices are assigned to match the index in the list
 * holding the molecule/atom. If a molecule is deleted during the course of the
 * simulation (e.g. grand-canonical sampling), the indices of other atoms are reassigned
 * to ensure they are contiguous; this circumstance triggers an atomLeafIndexChanged event,
 * to be acted upon by any registered listeners that may be affected by the reassignment.
 * <br>
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

    private static final long serialVersionUID = 2L;
    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList;
    protected final AtomSetAllMolecules allMoleculeList;
    private final BoxEventManager eventManager;
    private final ISpace space;
    protected int[] indexReservoir;
    protected int reservoirSize = 50;
    protected int reservoirCount;
    protected int maxIndex;
    protected MoleculeArrayList[] moleculeLists;
    private IBoundary boundary;
    private int index;
  
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
    
    public int getIndex() {
        return index;
    }

    public void setIndex(int newIndex) {
        index = newIndex;
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

    public final IBoundary getBoundary() {return boundary;}
    
    public void setBoundary(IBoundary b) {
        boundary = b;
        boundary.setBox(this);
     }

    public void setDensity(double rho) {
        double vNew = getMoleculeList().getMoleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        BoxInflate inflater = new BoxInflate(this, space);
        inflater.setScale(scale);
        inflater.actionPerformed();
    };

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
}
