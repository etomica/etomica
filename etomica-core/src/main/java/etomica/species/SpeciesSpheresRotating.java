/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeOriented;
import etomica.chem.elements.IElement;
import etomica.config.ConformationLinear;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOriented
 */
public class SpeciesSpheresRotating {

    public static SpeciesGeneral create(Space space, IElement element) {
        return create(space, element, false, true);
    }

    public static SpeciesGeneral create(Space space, IElement element, boolean isAxisSymmetric) {
        return create(space, element, false, isAxisSymmetric);
    }

    public static SpeciesGeneral create(Space space, IElement element, boolean isDynamic, boolean isAxisSymmetric) {
        Vector moment = space.makeVector();
        moment.E(1);
        if (isAxisSymmetric) {
            moment.setX(0, 0);
        }
        return new SpeciesBuilder(space)
                .setDynamic(isDynamic)
                .addCount(new AtomTypeOriented(element, moment), 1)
                .withConformation(new ConformationLinear(space))
                .withAtomFactory((atomType -> {
                    return isDynamic ? new AtomOrientedDynamic(space, atomType, isAxisSymmetric)
                            : new AtomOriented(space, atomType, isAxisSymmetric);
                }))
                .build();
    }
}
