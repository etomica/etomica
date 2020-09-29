/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.IConformation;
import etomica.meta.annotations.IgnoreProperty;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.List;

/**
 * An ISpecies holds information about how to construct a molecule, and creates molecules on demand.  In addition to
 * actually creating molecules, ISpecies objects are often used as keys to pass and designate which species is of
 * interest (for assigning potentials or setting the number of molecules in the box, for instance).
 *
 * @see Box
 * @see PotentialMaster
 */
public interface ISpecies {

    /**
     * Returns the index for this ISpecies, within the context of an Simulation.  The index is the ISpecies' position in
     * the array returned by SpeciesManager.getSpecies().
     */
    int getIndex();

    /**
     * Informs the ISpecies what its index should be.  This should only be called by the SpeciesManager.
     */
    void setIndex(int newIndex);

    IMolecule initMolecule(Box box, int molIdx, int atomIdxStart);

    /**
     * Returns whether the Species is dynamic. If a Species is dynamic and has oriented molecules, then the molecules
     * will have a velocity component. If a Species is dynamic and is not oriented, its atoms will have a velocity
     * component.
     *
     * @return whether the Species is dynamic
     */
    boolean isDynamic();

    /**
     *
     * @return the Space of the Species. This is the Space that molecules and atoms are constructed with.
     */
    Space getSpace();

    /**
     * Returns the number of unique AtomTypes that make up this Species.
     */
    @IgnoreProperty
    int getUniqueAtomTypeCount();

    /**
     *
     * @return the total number of atoms in this Species.
     */
    int getLeafAtomCount();

    /**
     * Returns the child type of this Species for the specified index.
     */
    AtomType getAtomType(int index);

    /**
     *
     * @return the singular AtomType that makes up this Species.
     */
    AtomType getLeafType();

    /**
     * Get the list of unique AtomTypes that make up this Species.
     * <p>
     * The returned list should not be modified.
     */
    List<AtomType> getUniqueAtomTypes();

    /**
     * Get a list with the AtomType for each atom which makes up this Species. The list will be in the same order as
     * the atoms within the molecule.
     * <p>
     * The returned list should not be modified.
     *
     * @return the list of AtomTypes of the atoms in this Species.
     */
    List<AtomType> getAtomTypes();

    /**
     * Initializes the molecule's IAtoms in their standard arrangement. The overall position of the molecule may be
     * changed by this method.
     */
    void initializeConformation(IMolecule molecule);

    /**
     *
     * @return the Conformation used by this Species to initialize positions of atoms within a molecule.
     */
    IConformation getConformation();

    /**
     * Get the index within a molecule instance of a named atom for this Species.
     *
     * @param atomName the name for a specific Atom in this Species.
     * @return an index for use in {@link IMolecule#getChildList()}
     */
    int getByName(String atomName);

    /**
     * Get the index within the molecule of the nth atom of given AtomType name.
     *
     * @param name String matching the name for the AtomType.
     * @param number The nth atom with that type to get, starting at 1.
     * @return the index of that atom within a molecule of this species.
     */
    int getAtomByTypeName(String name, int number);

    /**
     * Get the index within the molecule of the first atom of a given AtomType name.
     * @param name String matching the name for the AtomType.
     * @return the index of that atom within a molecule of this species.
     */
    int getAtomByTypeName(String name);

    /**
     * Get the AtomType with a given name used in this Species.
     * @param typeName String matching the name for the AtomType.
     * @return AtomType with the given name.
     */
    AtomType getTypeByName(String typeName);

    /**
     *
     * @return the moment of inertia for a molecule of this Species.
     */
    Vector getMomentOfInertia();

    /**
     *
     * @return the total mass of all the atoms in a molecule of this Species.
     */
    double getMass();
}
