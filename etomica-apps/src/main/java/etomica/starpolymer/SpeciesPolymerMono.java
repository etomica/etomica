/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.atom.AtomType;
import etomica.space.Space;
import etomica.species.SpeciesBuilder;

/**
 * Species in which molecules are each made of a single spherical atom.
 * Does not permit multiatomic molecules.  The advantage of this species
 * over the multiatomic version (used with 1 atom), is that one layer of
 * the atom hierarchy is eliminated in SpeciesSpheresMono.  Each atom is
 * the direct child of the species agent (i.e., each atom is at the "molecule"
 * level in the hierarchy, without an intervening AtomGroup).
 *
 * @author David Kofke
 */
public class SpeciesPolymerMono {

    public static SpeciesBuilder create(Space space, AtomType type, int f, int l) {
        return new SpeciesBuilder(space)
                .addCount(type, f * l + 1)
                .withConformation(new ConformationStarPolymerGraft(space, f, l));
    }
}
