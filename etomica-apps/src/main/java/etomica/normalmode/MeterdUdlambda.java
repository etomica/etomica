/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Energy;

public class MeterdUdlambda extends DataSourceScalar {

    protected final PotentialCompute pmModel, pmField;
    protected Box box;

    public MeterdUdlambda(PotentialCompute pmModel, PotentialCompute pmField, Box box) {
        super("dU/dlambda", Energy.DIMENSION);
        this.pmModel = pmModel;
        this.pmField = pmField;
        this.box = box;
    }

    public double getDataAsScalar() {
        pmModel.init();
//        pmField.init();
        double uModel = pmModel.computeAll(false);
//        double uField = pmField.computeAll(false);

//        double vx = ((IAtomKinetic) box.getLeafList().get(0)).getVelocity().getX(0);
//        double fx = pmModel.getForces()[0].getX(0);
//        System.out.println(uModel);

//        return uModel - uField;
//        System.out.println(uModel/box.getLeafList().size());
        return uModel;
    }

}