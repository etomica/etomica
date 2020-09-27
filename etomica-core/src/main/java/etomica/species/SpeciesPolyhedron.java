/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomTypeOriented;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.config.ConformationLinear;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.List;

/**
 * Species in which molecules are made of a single atom of type SpheroPolyhedron
 *
 * @author Andrew Schultz
 * @see AtomTypeOriented
 */
public class SpeciesPolyhedron {

    public static SpeciesGeneral create(Space space, List<Vector> vertices, double sweepRadius, IElement element) {
        return new SpeciesBuilder(space)
                .addCount(new AtomTypeSpheroPolyhedron(element, space, vertices, sweepRadius), 1)
                .withConformation(new ConformationLinear(space))
                .withAtomFactory(((atomType, box, id) -> {
                    return new AtomOrientedQuaternion(space, atomType);
                }))
                .build();
    }
}
