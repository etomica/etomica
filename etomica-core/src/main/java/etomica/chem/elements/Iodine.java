/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Chlorine element. 
 *
 * @author Andrew
 */
public class Iodine extends ElementChemical {

	protected Iodine(String symbol) {
        this(symbol, 126.9044);
    }

    protected Iodine(String symbol, double mass) {
        super(symbol, mass, 53);
    }

    public static final Iodine INSTANCE = new Iodine("I");
    
	private static final long serialVersionUID = 1L;
}
