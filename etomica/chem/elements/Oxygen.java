/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Oxygen element. 
 *
 * @author zhaofang
 */
public class Oxygen extends ElementChemical {
    
    protected Oxygen(String symbol) {
        this(symbol, 15.9994);
    }
    
    protected Oxygen(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    private static final long serialVersionUID = 1L;
    public static final Oxygen INSTANCE = new Oxygen("O");
}
