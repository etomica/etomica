/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.simulation.prototypes;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Angle;

public class MeterTorsionAngle extends DataSourceScalar {
    protected final Box box;
    protected final int[] atomsIdx;
    public MeterTorsionAngle(Box box, int iAtom, int jAtom, int kAtom, int lAtom) {
        super("torsion angle", Angle.DIMENSION);
        this.box = box;
        atomsIdx = new int[]{iAtom, jAtom, kAtom, lAtom};
    }

    @Override
    public double getDataAsScalar() {
        IAtomList atoms = box.getLeafList();
        Vector ri = atoms.get(atomsIdx[0]).getPosition();
        Vector rj = atoms.get(atomsIdx[1]).getPosition();
        Vector rk = atoms.get(atomsIdx[2]).getPosition();
        Vector rl = atoms.get(atomsIdx[3]).getPosition();
        Vector rji = new Vector3D();
        Vector rjk = new Vector3D();
        Vector rkl = new Vector3D();

        rji.Ev1Mv2(ri, rj);
        box.getBoundary().nearestImage(rji);
        rjk.Ev1Mv2(rk, rj);
        box.getBoundary().nearestImage(rjk);
        double rjk2 = rjk.squared();
        rkl.Ev1Mv2(rl, rk);
        box.getBoundary().nearestImage(rkl);

        Vector vji = new Vector3D();
        vji.E(rji);
        vji.PEa1Tv1(-rjk.dot(rji) / rjk2, rjk);
        double vji2 = vji.squared();
        Vector vkl = new Vector3D();
        vkl.E(rkl);
        vkl.PEa1Tv1(-rjk.dot(rkl) / rjk2, rjk);
        double vkl2 = vkl.squared();
        double rji2 = rji.squared();
        double rkl2 = rkl.squared();

        double vji2vkl2 = vji2 * vkl2;
        if (vji2 < 1e-6 * rji2 || vkl2 < 1e-6 * rkl2) {
            // one of the vectors (ji, kl) is nearly colinear with jk
            return 0;
        }
        double vji_vkl = 1 / Math.sqrt(vji2vkl2);
        return vji.dot(vkl) * vji_vkl;
    }
}
