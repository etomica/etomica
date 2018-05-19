/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.meta.annotations.IgnoreProperty;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;

import java.util.List;

/**
 * An ISpecies holds information about how to construct a molecule, and creates molecules on demand.  In addition to
 * actually creating molecules, ISpecies objects are often used as keys to pass and designate which species is of
 * interest (for assigning potentials or setting the number of molecules in the box, for instance).
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 * @author Andrew Schultz
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

    /**
     * Builds and returns the IMolecule of this ISpecies.
     */
    IMolecule makeMolecule();

    /**
     * Returns the number of child types of this group.
     */
    @IgnoreProperty
    int getAtomTypeCount();

    /**
     * Returns the child types of this group for the specified index.
     */
    AtomType getAtomType(int index);

    /**
     * Get the list of AtomTypes that make up this Species.
     * <p>
     * In general this list will be unmodifiable (will throw an exception if you attempt to modify it), instead you
     * should use methods defined on the implementing class to add or remove AtomTypes.
     *
     * @return
     */
    List<AtomType> getAtomTypes();

    /**
     * Initializes the molecule's IAtoms in their standard arrangement. The overall position of the molecule may be
     * changed by this method.
     */
    void initializeConformation(IMolecule molecule);
}
