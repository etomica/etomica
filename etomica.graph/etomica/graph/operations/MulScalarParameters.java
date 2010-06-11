package etomica.graph.operations;

import etomica.graph.model.Coefficient;

public class MulScalarParameters implements Parameters {

    private final Coefficient factor;

    public MulScalarParameters(Coefficient factor) {
        this.factor = factor;
    }

    public Coefficient factor() {

      return factor;
    }

}
