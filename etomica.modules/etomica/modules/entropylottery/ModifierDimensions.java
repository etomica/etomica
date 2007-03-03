package etomica.modules.entropylottery;

import etomica.modifier.Modifier;
import etomica.phase.Phase;
import etomica.space.IVectorRandom;
import etomica.units.Dimension;
import etomica.units.Length;

public class ModifierDimensions implements Modifier {

    public ModifierDimensions(Phase phase) {
        this.phase = phase;
    }

    public Dimension getDimension() {
        return Length.DIMENSION;
    }

    public String getLabel() {
        return "Dimensions";
    }

    public double getValue() {
        return phase.getBoundary().getDimensions().x(0);
    }

    public void setValue(double newValue) {
        if (newValue <= 0 || newValue > 1000) {
            throw new IllegalArgumentException("Bogus value for dimension");
        }
        IVectorRandom dim = phase.getBoundary().getDimensions();
        dim.setX(0, newValue);
        phase.setDimensions(dim);
    }

    private final Phase phase;
}
