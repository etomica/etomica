/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Helium element.
 */
public class Helium extends ElementChemical {
	
    protected Helium(String symbol) {
        this(symbol, 4.002602);
    }
	
    protected Helium(String symbol, double mass) {
        super(symbol, mass, 2);
    }
    
    public static final Helium INSTANCE = new Helium("He");
}
