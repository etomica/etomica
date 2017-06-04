/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.potential.PotentialMaster;

/**
 * An ISpecies holds information about how to construct a molecule, and
 * creates molecules on demand.  In addition to actually creating molecules,
 * ISpecies objects are often used as keys to pass and designate which
 * species is of interest (for assigning potentials or setting the number of
 * molecules in the box, for instance).
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 * @author Andrew Schultz
 * @see Box
 * @see PotentialMaster
 */
public interface ISpecies {

    /**
     * Informs the ISpecies what its index should be.  This should only be
     * called by the SpeciesManager.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the index for this ISpecies, within the context of an
     * ISimulation.  The index is the ISpecies' position in the array returned
     * by SpeciesManager.getSpecies().
     */
    public int getIndex();

    /**
     * Builds and returns the IMolecule of this ISpecies.
     */
    public IMolecule makeMolecule();

    /**
     * Returns the number of child types of this group.
     */
    public int getAtomTypeCount();

    /**
     * Returns the child types of this group for the specified index.
     */
    public IAtomType getAtomType(int index);

    /**
     * Initializes the molecule's IAtoms in their standard arrangement.
     * The overall position of the molecule may be changed by this method.
     */ 
    public void initializeConformation(IMolecule molecule);
}
