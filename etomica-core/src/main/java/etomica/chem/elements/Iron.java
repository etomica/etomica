/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Iron element.
 * reference: Atomic weights of the elements 2009 (IUPAC Technical Report)
 * http://iupac.org/publications/pac/83/2/0359/
 * 
 * @author Andrew Schultz
 */
public class Iron extends ElementChemical {

    protected Iron(String symbol) {
        this(symbol, 55.845);
    }
    
    public Iron(String symbol, double mass) {
        super(symbol, mass, 26);
    }

    public static final Iron INSTANCE = new Iron("Fe");
}
