/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Silver element.
 *
 * @author Andrew
 */
public class Silver extends ElementChemical {
    
    protected Silver(String symbol) {
        this(symbol, 107.8682);
    }
    
    public Silver(String symbol, double mass) {
        super(symbol, mass, 47);
    }

    private static final long serialVersionUID = 1L;
    public static final Silver INSTANCE = new Silver("Ag");
}
