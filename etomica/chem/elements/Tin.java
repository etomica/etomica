/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Tin element.
 *
 * @author Andrew
 */
public class Tin extends ElementChemical {
    
    protected Tin(String symbol) {
        this(symbol, 118.710);
    }
    
    public Tin(String symbol, double mass) {
        super(symbol, mass, 50);
    }

    private static final long serialVersionUID = 1L;
    public static final Tin INSTANCE = new Tin("Sn");
}
