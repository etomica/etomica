/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

/**
 * This class hashes atomic diameters based on both element and atom type.
 * Specifying diameters by type overrides any element specification (since type
 * can be more specific).
 *
 * @author Andrew Schultz
 */
public class DiameterHashByElementType extends DiameterHashByType {

    private final DiameterHashByElement diameterManagerByElement;

    public DiameterHashByElementType() {
        super();
        diameterManagerByElement = new DiameterHashByElement();
    }

    public double getDiameter(IAtom atom) {
        double d = super.getDiameter(atom);
        if (d >= 0) {
            return d;
        }
        return diameterManagerByElement.getDiameter(atom);
    }

    public void setDiameter(String element, double newDiameter) {
        diameterManagerByElement.setDiameter(element, newDiameter);
    }

    public DiameterHashByElement getDiameterHashByElement() {
        return diameterManagerByElement;
    }
}
