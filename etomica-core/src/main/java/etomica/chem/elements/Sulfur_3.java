/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Sulfur element.
 *
 * @author Andrew
 */
public class Sulfur_3 extends ElementChemical {

    protected Sulfur_3(String symbol) {
        this(symbol, 32.065);
    }

    public Sulfur_3(String symbol, double mass) {
        super(symbol, mass, 16);
    }

    private static final long serialVersionUID = 1L;
    public static final Sulfur_3 INSTANCE = new Sulfur_3("S_3");
}
