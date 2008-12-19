package etomica.box;

import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class RandomPositionSourceRectangular implements RandomPositionSource {
    
    public RandomPositionSourceRectangular(ISpace space, IRandom random) {
        p = (IVectorRandom)space.makeVector();
        this.random = random;
    }

    public IVectorMutable randomPosition() {
        p.setRandomCube(random);
        p.TE(box.getBoundary().getDimensions());
        return p;
    }

    public void setBox(IBox newBox) {
        box = newBox;
    }

    protected final IRandom random;
    protected IBox box;
    protected final IVectorRandom p;
}
