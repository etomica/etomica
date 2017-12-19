/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.util.HashMap;
import java.util.Map;

/**
 * This class hashes atomic diameters based on the atom type.
 *
 * @author Andrew Schultz
 */
public class DiameterHashByType implements DiameterHash {

    private final Map<AtomType, Double> diameterMap;

    public DiameterHashByType() {
        diameterMap = new HashMap<>();
    }

    public double getDiameter(IAtom atom) {
        return getDiameter(atom.getType());
    }

    public double getDiameter(AtomType atomType) {
        return diameterMap.getOrDefault(atomType, -1.0);
    }

    public void setDiameter(AtomType type, double newDiameter) {
        diameterMap.put(type, newDiameter);
    }
}
