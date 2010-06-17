package etomica.graph.operations;

import etomica.graph.model.Coefficient;
import etomica.graph.model.impl.CoefficientImpl;

public class MulScalarParameters implements Parameters {

    private final Coefficient factor;

    public MulScalarParameters(int numerator, int denominator) {
      factor = new CoefficientImpl(numerator, denominator);
    }

    public MulScalarParameters(Coefficient factor) {
        this.factor = factor;
    }

    public Coefficient factor() {

      return factor;
    }

}
