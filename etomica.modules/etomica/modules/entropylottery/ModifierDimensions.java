package etomica.modules.entropylottery;

import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.modifier.Modifier;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Length;

public class ModifierDimensions implements Modifier {

    public ModifierDimensions(ISpace space, IBox box) {
        this.box = box;
        boxDim = space.makeVector();
    }

    public Dimension getDimension() {
        return Length.DIMENSION;
    }

    public String getLabel() {
        return "Dimensions";
    }

    public double getValue() {
        return box.getBoundary().getBoxSize().getX(0);
    }

    public void setValue(double newValue) {
        if (newValue <= 0 || newValue > 1000) {
            throw new IllegalArgumentException("Bogus value for dimension");
        }
        boxDim.E(box.getBoundary().getBoxSize());
        boxDim.setX(0, newValue);
        box.getBoundary().setBoxSize(boxDim);
    }

    private final IBox box;
    protected final IVectorMutable boxDim;
}
