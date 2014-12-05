/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Chlorine element. 
 *
 * @author Andrew
 */
public class Chlorine extends ElementChemical {

	protected Chlorine(String symbol) {
        this(symbol, 35.453);
    }
    
    protected Chlorine(String symbol, double mass) {
        super(symbol, mass, 17);
    }

    public static final Chlorine INSTANCE = new Chlorine("Cl");
    
	private static final long serialVersionUID = 1L;
}
