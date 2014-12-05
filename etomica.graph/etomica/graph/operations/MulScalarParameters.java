/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
