/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Tungsten element.
 *
 * @author Andrew
 */
public class Tungsten extends ElementChemical {
    
    protected Tungsten(String symbol) {
        this(symbol, 183.84);
    }
    
    public Tungsten(String symbol, double mass) {
        super(symbol, mass, 74);
    }

    private static final long serialVersionUID = 1L;
    public static final Tungsten INSTANCE = new Tungsten("W");
}
