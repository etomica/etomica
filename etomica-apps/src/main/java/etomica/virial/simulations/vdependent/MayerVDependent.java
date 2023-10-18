/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.vdependent;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.virial.MayerFunction;

/**
 * This class wraps a spherical MayerFunction.  The value returned here is the
 * wrapped value after accounting for periodic boundaries.
 */
public class MayerVDependent implements MayerFunction {

    protected final Boundary boundary;
    protected final Vector dr;
    protected final MayerFunction f;

    protected MayerVDependent(MayerFunction f, Boundary b) {
        this.f = f;
        boundary = b;
        dr = boundary.getBox().getSpace().makeVector();
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        dr.Ev1Mv2(pair.get(0).getChildList().get(0).getPosition(),
                pair.get(1).getChildList().get(0).getPosition());
        boundary.nearestImage(dr);
        r2 = dr.squared();
        return f.f(pair, r2, beta);
    }

    public void setBox(Box newBox) {
    }
}
