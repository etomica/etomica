/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Mass;

/**
 * Abstract structure for a class defining an element.
 */
public abstract class Element implements IElement, java.io.Serializable {
	
	public Element(String symbol) {
        this.symbol = symbol;
	}
    
    public abstract double getMass();
    
    public Dimension getMassDimension() {
        return Mass.DIMENSION;
    }
    
    public abstract double rm();

    public String getSymbol() {
        return symbol;
    }
    
    protected final String symbol;
}
