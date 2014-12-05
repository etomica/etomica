/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Nitrogen element. 
 *
 * @author Andrew
 */
public class Nitrogen extends ElementChemical {
  
	protected Nitrogen(String symbol) {
        this(symbol, 14.01);
    }
    
    protected Nitrogen(String symbol, double mass) {
        super(symbol, mass, 7);
    }

    public static final Nitrogen INSTANCE = new Nitrogen("N");
    
	private static final long serialVersionUID = 1L;
}
