/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.UnitSystem;

/**
 * Class to form a dimension from ratio of two dimensions. Used primarily to
 * construct intensive quantities, such as energy/volume.
 */
public final class DimensionRatio extends Dimension {

    public DimensionRatio(Dimension nDimension, Dimension dDimension) {
        this(nDimension.toString() + "/" + dDimension.toString(), nDimension, dDimension);
    }

    public DimensionRatio(String name, Dimension nDimension, Dimension dDimension) {
        super(name, makeSignature(nDimension, dDimension));
        this.nDimension = nDimension;
        this.dDimension = dDimension;
    }

    private static double[] makeSignature(Dimension nDim, Dimension dDim) {
        double[] sig = nDim.signature().clone();
        for (int i = 0; i < sig.length; i++) {
            sig[i] -= dDim.signature()[i];
        }
        return sig;
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return new UnitRatio(nDimension.getUnit(unitSystem), dDimension.getUnit(unitSystem));
    }

    public Dimension nDimension() {
        return nDimension;
    }

    public Dimension dDimension() {
        return dDimension;
    }

    private final Dimension nDimension;
    private final Dimension dDimension;
    private static final long serialVersionUID = 1;

}
