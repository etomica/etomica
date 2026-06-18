/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterPressureHardTensor implements IDataSource, IntegratorHard.CollisionListener {

    private double t0;
    protected final Tensor v;
    protected final IntegratorHard integratorHard;
    private final Tensor virialSum;
    private final DataTensor data;
    private final IDataInfo dataInfo;
    protected final DataTag tag;
    protected final boolean justVirial;
    protected boolean doNonEquilibrium;
    protected boolean returned;

    public MeterPressureHardTensor(IntegratorHard integrator) {
        this(integrator, false);
    }

    public MeterPressureHardTensor(IntegratorHard integrator, boolean justVirial) {
        this.integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        this.justVirial = justVirial;
        Space space = integrator.getBox().getSpace();
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("pressure tensor", Pressure.dimension(space.D()), space);
        v = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
        virialSum = space.makeTensor();
        doNonEquilibrium = true;
        returned = true;
    }

    public void setDoNonEquilibrium(boolean doNonEquilibrium) {
        this.doNonEquilibrium = doNonEquilibrium;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        double t = integratorHard.getCurrentTime();
        if (returned) {
            // once we return a value, we need to zero out for new data
            data.x.E(0);
        }
        data.x.PEa1Tt1(-1 / ((t - t0)), virialSum);
        virialSum.E(0);
        t0 = t;

        if (justVirial) {
            returned = true;
            return data;
        }

        //We're using the instantaneous velocity tensor with the average virial tensor
        //not quite right, but (so long as you need averages) it doesn't really matter
        // if you want more "correct" properties (have both be an average), then you'll
        // need to compute a correction in collisionAction based on the old and new
        // velocities
        Box box = integratorHard.getBox();
        IAtomList leafList = box.getLeafList();
        if (doNonEquilibrium) {
            // We're using the current velocity tensor with the average virial tensor
            // We've computed corrections to this within pairCollision
            Tensor K = box.getSpace().makeTensor();
            for (IAtom iAtom : leafList) {
                Vector v = ((IAtomKinetic) iAtom).getVelocity();
                K.Ev1v2(v, v);
                K.TE(iAtom.getType().getMass());
                data.x.PE(K);
            }
        }
        else {
            // include an average value for the kinetic contribution
            Vector v = box.getSpace().makeVector();
            v.E(1);
            Tensor I = box.getSpace().makeTensor();
            I.diagE(v);
            data.x.PEa1Tt1(leafList.size() * integratorHard.getTemperature(), I);
        }

        data.x.TE(1.0 / box.getBoundary().volume());
        returned = true;

        return data;
    }

    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector dv, double virial, double tCollision) {
        v.Ev1v2(r12, r12);
        v.TE(virial / r12.squared());
        virialSum.PE(v);

        if (returned) {
            // we're starting fresh
            data.x.E(0);
            returned = false;
        }

        if (doNonEquilibrium && !justVirial) {
            // In getData, we'll estimate the kinetic contributions based on the velocities at the
            // end of the step. Right now, we need to compute a correction to that estimate based
            // on how much the velocities changed due to this collision.
            // Basically, since the beginning of the timestep, we have had the pre-collision
            // velocities.  So include their contribution, but also subtract off the contribution
            // from the new velocities (which we will incorrectly include later).
            // If an atom happens to have multiple collisions in a step, all this stuff will
            // still work out.

            // FIXME This does depend on getData being called every step because our tCollision here
            // FIXME is only the time since the start of the current step.
            Space space = integratorHard.getBox().getSpace();
            double lastCollisionVirialr2 = virial / r12.squared();
            Vector dp = space.makeVector();
            dp.Ea1Tv1(lastCollisionVirialr2, r12);

            Vector dva = space.makeVector();
            dva.Ea1Tv1(1/atom1.getType().rm(), dp);
            Vector vOld = space.makeVector();
            vOld.Ev1Mv2(atom1.getVelocity(), dva);

            Tensor tTmp = space.makeTensor();
            Tensor tCor = space.makeTensor();
            tTmp.Ev1v2(vOld, vOld);
            tCor.PEa1Tt1(atom1.getType().getMass(), tTmp);
            tTmp.Ev1v2(atom1.getVelocity(), atom1.getVelocity());
            tCor.PEa1Tt1(-atom1.getType().getMass(), tTmp);

            dva.Ea1Tv1(-1/atom2.getType().rm(), dp);
            vOld.Ev1Mv2(atom2.getVelocity(), dva);

            tTmp.Ev1v2(vOld, vOld);
            tCor.PEa1Tt1(atom2.getType().getMass(), tTmp);
            tTmp.Ev1v2(atom2.getVelocity(), atom2.getVelocity());
            tCor.PEa1Tt1(-atom2.getType().getMass(), tTmp);

            data.x.PEa1Tt1(tCollision, tCor);
        }
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
}
