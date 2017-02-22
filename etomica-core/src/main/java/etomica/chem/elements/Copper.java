/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Copper element.
 *
 * @author Andrew
 */
public class Copper extends ElementChemical {
    
    protected Copper(String symbol) {
        this(symbol, 63.546);
    }
    
    public Copper(String symbol, double mass) {
        super(symbol, mass, 29);
    }

    private static final long serialVersionUID = 1L;
    public static final Copper INSTANCE = new Copper("Cu");
}
