package etomica.chem.elements;

import etomica.atom.AtomTypeRoot;
import etomica.simulation.Simulation;

public class ElementSimple extends Element {

    public ElementSimple(Simulation sim) {
        this(((AtomTypeRoot)sim.speciesRoot.getType()).makeUniqueElementSymbol("E"), sim.getDefaults().atomMass);
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
