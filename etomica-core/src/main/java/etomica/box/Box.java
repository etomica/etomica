/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.action.BoxInflate;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.storage.*;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Arrays;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
    protected final MoleculeArrayList allMoleculeList;
    private final BoxEventManager eventManager;
    private final Space space;
    protected MoleculeArrayList[] moleculeLists;
    private final Boundary boundary;
    private int index;

    private int atomCount;
    private int moleculeCount;
    private final Map<Object, VectorStorage> atomVectorStorageMap = new HashMap<>();
    private final Map<Object, DoubleStorage> atomDoubleStorageMap = new HashMap<>();
    private final Map<Object, OrientationStorage> atomOrientationsStorageMap = new HashMap<>();
    private final Map<Object, VectorStorage> molVectorStorageMap = new HashMap<>();
    private final Map<Object, DoubleStorage> molDoubleStorageMap = new HashMap<>();
    private final Map<Object, OrientationStorage> molOrientationsStorageMap = new HashMap<>();

    private final List<Storage> allAtomStorage = new ArrayList<>();
    private final List<Storage> allMoleculeStorage = new ArrayList<>();


    public VectorStorage getAtomVectors(Object token) {
        return this.atomVectorStorageMap.computeIfAbsent(token, t -> {
            VectorStorage storage = new VectorStorage(space, atomCount);
            allAtomStorage.add(storage);
            return storage;
        });
    }

    public DoubleStorage getAtomDoubles(Object token) {
        return this.atomDoubleStorageMap.computeIfAbsent(token, t -> {
            DoubleStorage storage = new DoubleStorage(atomCount);
            allAtomStorage.add(storage);
            return storage;
        });
    }

    public OrientationStorage getAtomOrientations(Tokens.OrientationsToken token) {
        boolean isAxisSymmetric = token.isAxisSymmetric();
        return this.atomOrientationsStorageMap.computeIfAbsent(token, t -> {
            OrientationStorage storage = new OrientationStorage(space, atomCount, isAxisSymmetric);
            allAtomStorage.add(storage);
            return storage;
        });
    }

    public VectorStorage getMolVectors(Object token) {
        return this.molVectorStorageMap.computeIfAbsent(token, t -> {
            VectorStorage storage = new VectorStorage(space, moleculeCount);
            allMoleculeStorage.add(storage);
            return storage;
        });
    }

    public DoubleStorage getMolDoubles(Object token) {
        return this.molDoubleStorageMap.computeIfAbsent(token, t -> {
            DoubleStorage storage = new DoubleStorage(moleculeCount);
            allMoleculeStorage.add(storage);
            return storage;
        });
    }

    public OrientationStorage getMolOrientations(Tokens.OrientationsToken token) {
        boolean isAxisSymmetric = token.isAxisSymmetric();
        return this.molOrientationsStorageMap.computeIfAbsent(token, t -> {
            OrientationStorage storage = new OrientationStorage(space, moleculeCount, isAxisSymmetric);
            allMoleculeStorage.add(storage);
            return storage;
        });
    }

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
        allMoleculeList = new MoleculeArrayList();
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
        return addNewMolecule(species, mol -> {});
    }

    public IMolecule addNewMolecule(ISpecies species, Consumer<IMolecule> initMolecule) {
        int numAtoms = species.getLeafAtomCount();

        int molIdx = this.addMoleculeStorage(1);
        int atomIdxStart = this.addAtomStorage(numAtoms);

        IMolecule molecule = species.initMolecule(this, molIdx, atomIdxStart);
        int speciesIndex = species.getIndex();
        molecule.setIndex(moleculeLists[speciesIndex].size());
        moleculeLists[speciesIndex].add(molecule);
        allMoleculeList.add(molecule);

        int nListAtoms = leafList.size();
        for (int i = 0; i < molecule.getChildList().size(); i++) {
            IAtom childAtom = molecule.getChildList().get(i);
            childAtom.setLeafIndex(nListAtoms++);
            leafList.add(childAtom);
        }

        initMolecule.accept(molecule);
        eventManager.moleculeAdded(molecule);
        return molecule;
    }

    public void removeMolecule(IMolecule molecule) {
        int molSpeciesIdx = molecule.getIndex();
        int molIdx = molecule.getGlobalIndex();
        MoleculeArrayList molList = moleculeLists[molecule.getType().getIndex()];
        if (molList.get(molSpeciesIdx) != molecule) {
            throw new IllegalArgumentException("can't find " + molecule);
        }

        molList.removeAndReplace(molSpeciesIdx);
        if (molSpeciesIdx < molList.size()) {
            IMolecule replacingMolecule = molList.get(molSpeciesIdx);
            replacingMolecule.setIndex(molSpeciesIdx);
            eventManager.moleculeIndexChanged(replacingMolecule, molList.size());
        }

        this.allMoleculeStorage.forEach(s -> s.swapRemove(molIdx));

        allMoleculeList.removeAndReplace(molIdx);
        if (molIdx < allMoleculeList.size()) {
            IMolecule replacingMolecule = allMoleculeList.get(molIdx);
            replacingMolecule.setGlobalIndex(molIdx);
        }

        eventManager.moleculeRemoved(molecule);

        for (IAtom atom : molecule.getChildList()) {
            removeAtom(atom.getLeafIndex());
        }

    }

    private void removeAtom(int atomIdx) {
        this.allAtomStorage.forEach(s -> s.swapRemove(atomIdx));
        this.atomCount--;
        leafList.removeAndReplace(atomIdx);
        if (atomIdx < leafList.size()) {
            IAtom replacingAtom = leafList.get(atomIdx);
            replacingAtom.setLeafIndex(atomIdx);
            eventManager.atomLeafIndexChanged(replacingAtom, leafList.size());
        }
    }

    public int addMoleculeStorage(int nMolecules) {
        this.allMoleculeStorage.forEach(s -> s.addNull(nMolecules));
        int newIdxStart = this.moleculeCount;
        this.moleculeCount += nMolecules;
        return newIdxStart;
    }

    public int addAtomStorage(int nAtoms) {
        this.allAtomStorage.forEach(s -> s.addNull(nAtoms));
        int newIdxStart = this.atomCount;
        this.atomCount += nAtoms;
        return newIdxStart;
    }

    private void ensureCapacityMolecules(int newMolecules) {
        this.allMoleculeStorage.forEach(s -> s.ensureCapacity(newMolecules));
    }

    private void ensureCapacityAtoms(int newAtoms) {
        this.allAtomStorage.forEach(s -> s.ensureCapacity(newAtoms));
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
        if (n < 0) {
            throw new IllegalArgumentException("Number of molecules cannot be negative");
        }
        int speciesIndex = species.getIndex();
        MoleculeArrayList moleculeList = moleculeLists[speciesIndex];
        int currentNMolecules = moleculeList.size();
        if (n > currentNMolecules) {
            int moleculeLeafAtoms = species.getLeafAtomCount();
            notifyNewMolecules(species, (n - currentNMolecules), moleculeLeafAtoms);
            moleculeLists[species.getIndex()].ensureCapacity(n);
            leafList.ensureCapacity(leafList.size() + (n - currentNMolecules) * moleculeLeafAtoms);

            ensureCapacityAtoms((n - currentNMolecules) * moleculeLeafAtoms);
            ensureCapacityMolecules(n - currentNMolecules);

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
        moleculeLists = Arrays.addObject(moleculeLists, new MoleculeArrayList());
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
