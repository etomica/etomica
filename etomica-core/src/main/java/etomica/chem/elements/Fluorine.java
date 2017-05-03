/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

/**
 * Class for the Fluorine element.
 * reference: Atomic weights of the elements 2009 (IUPAC Technical Report)
 * http://iupac.org/publications/pac/83/2/0359/
 * 
 * @author Shu
 */
public class Fluorine extends ElementChemical {

    protected Fluorine(String symbol) {
        this(symbol, 18.9984);
    }
    
    public Fluorine(String symbol, double mass) {
        super(symbol, mass, 9);
    }

    private static final long serialVersionUID = 1L;
    public static final Fluorine INSTANCE = new Fluorine("F");
}
