/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.atom.IAtomKinetic;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorHardFasterer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure2D;

public class DataSourceWallPressureMuFasterer implements IDataSource, IntegratorHardFasterer.CollisionListener {
    public DataSourceWallPressureMuFasterer(Space space) {
        this.space = space;
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("pressure", Pressure2D.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    @Override
    public void fieldCollision(IAtomKinetic atom, Vector r, Vector deltaP, double tCollision) {
        double x = r.getX(0);
        if (Math.abs(x) < 1) return;
        if (x < 0) {
            virialSumIG -= deltaP.getX(0);
        } else {
            virialSumSQW += deltaP.getX(0);
        }
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        double currentTime = integratorHard.getCurrentTime();
        double[] x = data.getData();
        x[0] = virialSumSQW / (currentTime - lastTime);
        x[1] = virialSumIG / (currentTime - lastTime);
        lastTime = currentTime;
        virialSumIG = 0;
        virialSumSQW = 0;
        return data;
    }

    /**
     * Registers meter as a collisionListener to the integrator, and sets up
     * a DataSourceTimer to keep track of elapsed time of integrator.
     */
    public void setIntegrator(IntegratorHardFasterer newIntegrator) {
        if (newIntegrator == integratorHard) return;
        if (integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
        if (newIntegrator != null) {
            integratorHard.addCollisionListener(this);
            lastTime = integratorHard.getCurrentTime();
        }
        virialSumIG = virialSumSQW = 0;
    }

    public IntegratorHardFasterer getIntegrator() {
        return integratorHard;
    }

    private static final long serialVersionUID = 1L;
    protected Space space;
    protected IntegratorHardFasterer integratorHard;
    protected double virialSumIG, virialSumSQW;
    protected double lastTime;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision) {

    }

}
