package etomica.modules.pistoncylinder;

import etomica.Modifier;
import etomica.potential.P1HardMovingBoundary;
import etomica.units.Dimension;

/**
 * Modifier for the piston pressure.  Should be followed by a IntegratorPistonUpdate action.
 */
class ModifierPistonPressure implements Modifier {
    public ModifierPistonPressure(P1HardMovingBoundary potential, Dimension pressureDimension) {
        pistonPotential = potential;
        pressureDim = pressureDimension;
    }

    public void setValue(double p) {
        pistonPotential.setPressure(p);
    }

    public double getValue() {
        return pistonPotential.getPressure();
    }

    public Dimension getDimension() {
        return pressureDim;
    }
    
    public String getLabel() {
        return "Piston pressure";
    }
    
    public String toString() {
        return getLabel();
    }
    private final P1HardMovingBoundary pistonPotential;
    private final Dimension pressureDim;
}