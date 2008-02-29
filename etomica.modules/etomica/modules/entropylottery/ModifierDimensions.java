package etomica.modules.entropylottery;

import etomica.modifier.Modifier;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.units.Dimension;
import etomica.units.Length;

public class ModifierDimensions implements Modifier {

    public ModifierDimensions(Box box) {
        this.box = box;
    }

    public Dimension getDimension() {
        return Length.DIMENSION;
    }

    public String getLabel() {
        return "Dimensions";
    }

    public double getValue() {
        return box.getBoundary().getDimensions().x(0);
    }

    public void setValue(double newValue) {
        if (newValue <= 0 || newValue > 1000) {
            throw new IllegalArgumentException("Bogus value for dimension");
        }
        IVector dim = box.getBoundary().getDimensions();
        dim.setX(0, newValue);
        box.setDimensions(dim);
    }

    private final Box box;
}
