/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Energy;

public class MeterdUdlambda extends DataSourceScalar {

    protected Box box;
    protected final PotentialCompute pmModel, pmField;

    public MeterdUdlambda(PotentialCompute pmModel, PotentialCompute pmField) {
        super("dU/dlambda", Energy.DIMENSION);
        this.pmModel = pmModel;
        this.pmField = pmField;
    }

    public double getDataAsScalar() {
        double uModel = pmModel.computeAll(false);
        double uField = pmField.computeAll(false);
        return uModel - uField;
    }

}
