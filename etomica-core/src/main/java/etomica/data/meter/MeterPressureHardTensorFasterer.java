/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHardFasterer;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterPressureHardTensorFasterer implements IDataSource, IntegratorHardFasterer.CollisionListener {

    public MeterPressureHardTensorFasterer(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("pressure tensor", Pressure.dimension(space.D()), space);
        v = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
        virialSum = space.makeTensor();
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        if (box == null || integratorHard == null)
            throw new IllegalStateException("must call setBox and integrator before using meter");
        double t = integratorHard.getCurrentTime();
        data.x.E(0);
        data.x.PEa1Tt1(-1 / ((t - t0)), virialSum);
        virialSum.E(0);
        t0 = t;

        //We're using the instantaneous velocity tensor with the average virial tensor
        //not quite right, but (so long as you need averages) it doesn't really matter
        // if you want more "correct" properties (have both be an average), then you'll
        // need to compute a correction in collisionAction based on the old and new
        // velocities
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            v.Ev1v2(a.getVelocity(), a.getVelocity());
            v.TE((a.getType().rm()));
            data.x.PE(v);
        }

        data.x.TE(1.0 / box.getBoundary().volume());

        return data;
    }

    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector dv, double virial, double tCollision) {
        v.Ev1v2(r12, r12);
        v.TE(virial / r12.squared());
        virialSum.PE(v);
    }

    public void setIntegrator(IntegratorHardFasterer newIntegrator) {
        if (newIntegrator == integratorHard) {
            return;
        }
        if (integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        if (newIntegrator == null) {
            box = null;
            return;
        }
        integratorHard = newIntegrator;
        box = integratorHard.getBox();
        integratorHard.addCollisionListener(this);
    }

    public IntegratorHardFasterer getIntegrator() {
        return integratorHard;
    }

    private double t0;
    private final Tensor v;
    private IntegratorHardFasterer integratorHard;
    private Box box;
    private final Tensor virialSum;
    private final DataTensor data;
    private final IDataInfo dataInfo;
    protected final DataTag tag;
}
