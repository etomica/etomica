/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Carbon element. 
 *
 * @author Andrew
 */
public class Carbon extends ElementChemical {

	protected Carbon(String symbol) {
        this(symbol, 12.0107);
    }
    
    protected Carbon(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon INSTANCE = new Carbon("C");
    
	private static final long serialVersionUID = 1L;
}
