/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.chem.elements;

import etomica.simulation.Simulation;

public class ElementSimple extends Element {

    public ElementSimple(Simulation sim) {
        this(sim.makeUniqueElementSymbol("E"));
    }
    
    public ElementSimple(String symbol) {
        this(symbol, 1.0);
    }
    
    public ElementSimple(String symbol, double mass) {
        super(symbol);
        setMass(mass);
    }
    
    /**
     * Sets mass of this element and updates reciprocal mass accordingly.
     */
    public void setMass(double newMass) {
        mass = newMass;
        rm = 1.0/mass;
    }
    
    public final double getMass() {
        return mass;
    }
    
    public final double rm() {
        return rm;
    }
    
    private static final long serialVersionUID = 1L;
    protected double rm;
    protected double mass;
}
