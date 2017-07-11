/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Quantity;

/**
 * Abstract structure for a class defining one of the chemical elements.
 * Subclasses are (or will be) defined to correspond to each element in the
 * periodic table.
 */
public class ElementChemical extends Element {

    public ElementChemical(String symbol, double mass, int atomicNumber) {
        super(symbol);
        this.atomicNumber = atomicNumber;
        this.mass = mass;
        rm = 1.0/mass;
    }
    
    public final double getMass() {
        return mass;
    }
    
    public final double rm() {
        return rm;
    }
    
    public int getAtomicNumber() {
        return atomicNumber;
    }
    
    public Dimension getAtomicNumberDimension() {
        return Quantity.DIMENSION;
    }

    protected final double mass;
    protected final double rm;
    protected final int atomicNumber;
}
