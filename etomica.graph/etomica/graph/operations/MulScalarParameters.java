package etomica.graph.operations;

import etomica.graph.model.Coefficient;
import etomica.graph.model.impl.CoefficientImpl;

public class MulScalarParameters implements Parameters {

    private final Coefficient factor;

    public MulScalarParameters(int numerator, int denominator) {
      int sign = 1;
      if (numerator*denominator < 0) {
        sign = -1;
      }
      if (numerator < 0) {
        numerator = -numerator;
      }
      if (denominator < 0) {
        denominator = -denominator;
      }
      factor = new CoefficientImpl(numerator, denominator, sign);
    }

    public MulScalarParameters(Coefficient factor) {
        this.factor = factor;
    }

    public Coefficient factor() {

      return factor;
    }

}
