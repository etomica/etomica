/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.action.BoxInflate;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Arrays;
import etomica.util.Debug;

import java.util.function.Consumer;

/**
 * A Box collects all atoms that interact with one another; atoms in different
 * boxes do not interact. These are the important features of a Box:
 * <p>
 * <ul>
 * <li>It holds lists of the molecules in the box, including a list of all
 * molecules, and an array of lists that can be selected by Species. It
 * contains methods to add and remove molecules.</li>
 * <li>It holds a list of all atoms in the box.</li>
 * <li>It holds a Boundary object that
 * defines the volume of the box and the behavior of atoms as they move into or
 * across the boundary of the box.</li>
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the box (such as when molecules are added or removed).</li>
 * <li>Each Box has a unique index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Box
 * construction. This index is useful when collecting things in reference to the
 * box.</li>
 * </ul>
 * The box maintains index values for each molecule and each atom. These values
 * are held as integer fields by the corresponding Molecule and Atom instances.
 * The indices are assigned when molecule (and its atoms) are added to the box via
 * the addMolecule method. The indices are assigned to match the index in the list
 * holding the molecule/atom.
 * <br>
 * A box is acted upon by an Integrator instance to move its atoms around and
 * generate configurations. <br>
 * A simulation may involve more than one box. All Box instances should be
 * registered with the simulation via {@link etomica.simulation.Simulation#addBox(Box)} and
 * may be accessed via the simulation's getBox method.
 *
 * @author David Kofke, Andrew Schultz
 * @see Boundary
 * @see BoxEventManager
 */
public class Box {

    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList;
    protected final AtomSetAllMolecules allMoleculeList;
    private final BoxEventManager eventManager;
    private final Space space;
    protected MoleculeArrayList[] moleculeLists;
    private final Boundary boundary;
    private int index;

    /**
     * Constructs box with default rectangular periodic boundary.
     *
     * @param space space governing the simulation
     */
    public Box(Space space) {
        this(new BoundaryRectangularPeriodic(space), space);
    }

    /**
     * Constructs box with the given Boundary, which specifies the size and shape of the Box and specifies what
     * happens to atoms as they cross the box boundary e.g. periodic boundaries.
     *
     * @param boundary the specified Boundary
     * @param space    space governing the simulation
     */
    public Box(Boundary boundary, Space space) {
        this.space = space;
        eventManager = new BoxEventManager(this);
        this.boundary = boundary;
        this.boundary.setBox(this);

        moleculeLists = new MoleculeArrayList[0];
        allMoleculeList = new AtomSetAllMolecules();
        allMoleculeList.setMoleculeLists(moleculeLists);
        leafList = new AtomArrayList();
    }

    public Space getSpace() {
        return space;
    }

    /**
     * @return the Box's index.  The index corresponds to the box's position
     * in the simulation's list of Boxes.  The index of the first Box is 0.
     * The index of the last Box is n-1, where n is the number of Boxes.
     */
    public int getIndex() {
        return index;
    }

    /**
     * Informs the Box what its index is.  This should only be called by the
     * Simulation.
     *
     * @param newIndex the Box's new index
     */
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    /**
     * @return the String "Box" with the box index appended
     */
    public String toString() {
        return "Box " + getIndex();
    }

    /**
     * Creates a molecule of the given species and adds it to the box
     *
     * @param species the given species
     * @return the new molecule
     */
    public IMolecule addNewMolecule(ISpecies species) {
        return this.addNewMolecule(species, (mol) -> {});
    }

    public IMolecule addNewMolecule(ISpecies species, Consumer<IMolecule> initMolecule) {
        IMolecule mol = species.makeMolecule();
        addMolecule(mol, initMolecule);
        return mol;
    }

    /**
     * Adds the given molecule to the this Box.  The molecule should not
     * already be in this Box and should not be in another Box.  The molecule
     * should be a member of an ISpecies which has been added to the
     * Simulation. No exceptions are thrown if these conditions are violated.
     *
     * @param molecule the molecule to be added to the Box
     */
    private void addMolecule(IMolecule molecule, Consumer<IMolecule> initMolecule) {
        int speciesIndex = molecule.getType().getIndex();
        if (Debug.ON) {
            for (int i = 0; i < moleculeLists[speciesIndex].size(); i++) {
                if (moleculeLists[speciesIndex].get(i) == molecule) {
                    throw new RuntimeException("you bastard!");
                }
            }
        }
        molecule.setIndex(moleculeLists[speciesIndex].size());
        moleculeLists[speciesIndex].add(molecule);
        allMoleculeList.setMoleculeLists(moleculeLists);

        IAtomList childList = molecule.getChildList();
        int nLeafAtoms = leafList.size();
        for (int iChild = 0; iChild < childList.size(); iChild++) {
            IAtom childAtom = childList.get(iChild);
            childAtom.setLeafIndex(nLeafAtoms++);
            leafList.add(childAtom);
        }
        initMolecule.accept(molecule);
        eventManager.moleculeAdded(molecule);

        if (Debug.ON) {
            for (int i = 0; i < moleculeLists[speciesIndex].size(); i++) {
                if (moleculeLists[speciesIndex].get(i).getIndex() != i) {
                    throw new RuntimeException("oops " + molecule + " " + moleculeLists[speciesIndex].get(i) + " " + i);
                }
            }
        }
    }

    /**
     * Removes the given molecule from this Box.
     *
     * @param molecule the molecule to be removed from the Box
     * @throws IllegalArgumentException if the given molecule is not in the Box
     */
    public void removeMolecule(IMolecule molecule) {
        int moleculeIndex = molecule.getIndex();
        MoleculeArrayList moleculeList = moleculeLists[molecule.getType().getIndex()];
        if (moleculeList.get(moleculeIndex) != molecule) {
            throw new IllegalArgumentException("can't find " + molecule);
        }
        if (moleculeIndex < moleculeList.size() - 1) {
            moleculeList.removeAndReplace(moleculeIndex);
            IMolecule replacingMolecule = moleculeList.get(moleculeIndex);
            replacingMolecule.setIndex(moleculeIndex);
            eventManager.moleculeIndexChanged(replacingMolecule, moleculeList.size());
        } else {
            moleculeList.remove(moleculeIndex);
        }
        allMoleculeList.setMoleculeLists(moleculeLists);

        eventManager.moleculeRemoved(molecule);
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild < childList.size(); iChild++) {
            IAtom childAtom = childList.get(iChild);
            int leafIndex = childAtom.getLeafIndex();
            leafList.removeAndReplace(leafIndex);
            if (leafList.size() > leafIndex) {
                IAtom movedAtom = leafList.get(leafIndex);
                int movedLeafIndex = movedAtom.getLeafIndex();
                movedAtom.setLeafIndex(leafIndex);
                eventManager.atomLeafIndexChanged(movedAtom, movedLeafIndex);
            }
        }
        leafList.maybeTrimToSize();
    }

    /**
     * Sets the number of molecules in this box of the given Species to n.
     * Molecules are added to or removed from the box to achieve the desired
     * number.
     *
     * @param species the species whose number of molecules should be changed
     * @param n       the desired number of molecules
     * @throws IllegalArgumentException if n < 0.
     */
    public void setNMolecules(ISpecies species, int n) {
        int speciesIndex = species.getIndex();
        MoleculeArrayList moleculeList = moleculeLists[speciesIndex];
        int currentNMolecules = moleculeList.size();
        int moleculeLeafAtoms = 0;
        if (currentNMolecules > 0) {
            moleculeLeafAtoms = moleculeList.get(0).getChildList().size();
        } else if (n > currentNMolecules) {
            moleculeLeafAtoms = species.getLeafAtomCount();
        }
        notifyNewMolecules(species, (n - currentNMolecules), moleculeLeafAtoms);
        if (n < 0) {
            throw new IllegalArgumentException("Number of molecules cannot be negative");
        }
        if (n > currentNMolecules) {
            moleculeLists[species.getIndex()].ensureCapacity(n);
            leafList.ensureCapacity(leafList.size() + (n - currentNMolecules) * moleculeLeafAtoms);
            for (int i = currentNMolecules; i < n; i++) {
                this.addNewMolecule(species);
            }
        } else {
            for (int i = currentNMolecules; i > n; i--) {
                removeMolecule(moleculeList.get(i - 1));
            }
        }
    }

    /**
     * @return the number of molecules in this box of the given ISpecies.
     */
    public int getNMolecules(ISpecies species) {
        int speciesIndex = species.getIndex();
        return moleculeLists[speciesIndex].size();
    }

    /**
     * Returns the list of molecules of the given species.
     *
     * @param species the species
     * @return the requested list of molecules
     */
    public IMoleculeList getMoleculeList(ISpecies species) {
        return moleculeLists[species.getIndex()];
    }

    /**
     * @return a list of all molecules in this box.
     */
    public IMoleculeList getMoleculeList() {
        return allMoleculeList;
    }

    /**
     * @return the box's boundary.
     */
    public final Boundary getBoundary() {
        return boundary;
    }

    /**
     * Uses BoxInflate to adjust the volume to the specified density.
     * New volume is set such that N/V = rho, where N is the number of
     * molecules in the box.
     *
     * @param rho the specified density
     */
    public void setDensity(double rho) {
        double vNew = getMoleculeList().size() / rho;
        double scale = Math.pow(vNew / boundary.volume(), 1.0 / space.D());
        BoxInflate inflater = new BoxInflate(this, space);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }

    /**
     * @return the event manager for this box.
     */
    public BoxEventManager getEventManager() {
        return eventManager;
    }

    /**
     * Notifies the Box that the given species has been added to the
     * simulation.  This method should only be called by the simulation.
     *
     * @param species the added species
     */
    public void addSpeciesNotify(ISpecies species) {
        moleculeLists = (MoleculeArrayList[]) Arrays.addObject(moleculeLists, new MoleculeArrayList());
        allMoleculeList.setMoleculeLists(moleculeLists);
    }

    /**
     * Notifies the Box that a Species has been removed.  This method should
     * only be called by the simulation.  This triggers the removal of all
     * molecules of the given species from this box.
     *
     * @param species the removed species
     */
    public void removeSpeciesNotify(ISpecies species) {
        moleculeLists = (MoleculeArrayList[]) Arrays.removeObject(moleculeLists, moleculeLists[species.getIndex()]);
        allMoleculeList.setMoleculeLists(moleculeLists);
    }

    /**
     * @return the list of atoms contained in this box.
     */
    public IAtomList getLeafList() {
        return leafList;
    }


    protected void notifyNewMolecules(ISpecies species, int numNewMolecules, int moleculeLeafAtoms) {
        if (numNewMolecules < 1) return;
        // has no actual effect within this object.  We just notify things to
        // prepare for an increase in the max index.  If things assume that the
        // actual max index has already increased, there's no harm since
        // there's nothing that says the max index can't be too large.
        int numNewLeafAtoms = numNewMolecules * moleculeLeafAtoms;
        eventManager.numberMolecules(species, moleculeLists[species.getIndex()].size() + numNewMolecules);
        if (numNewLeafAtoms > 1) {
            eventManager.globalAtomLeafIndexChanged(leafList.size() + numNewLeafAtoms);
        }
    }
}
