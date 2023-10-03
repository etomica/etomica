/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Hydrogen element. 
 *
 * @author zhaofang
 */
public class Hydrogen_p extends ElementChemical {

    protected Hydrogen_p(String symbol) {
        this(symbol, 1.00794);
    }

    protected Hydrogen_p(String symbol, double mass) {
        super(symbol, mass, 1);
    }
    
    private static final long serialVersionUID = 1L;
    public static final Hydrogen_p INSTANCE = new Hydrogen_p("H_p");
}
