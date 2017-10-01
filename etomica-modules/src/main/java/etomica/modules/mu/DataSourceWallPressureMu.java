/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorHard.CollisionListener;
import etomica.potential.P1HardBoundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure2D;

public class DataSourceWallPressureMu implements IDataSource, CollisionListener {
    public DataSourceWallPressureMu(Space space) {
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
    public void collisionAction(IntegratorHard.Agent agent) {
        if (agent.collisionPotential instanceof P1HardBoundary) {
            Vector p = agent.atom.getPosition();
            if (p.getX(0) < 0) {
                virialSumIG -= ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
            }
            else {
                virialSumSQW += ((P1HardBoundary)agent.collisionPotential).lastWallVirial();
            }
        }
    }

    public IEtomicaDataInfo getDataInfo() {
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
    public void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
        if(newIntegrator != null) {
            integratorHard.addCollisionListener(this);
            lastTime = integratorHard.getCurrentTime();
        }
        virialSumIG = virialSumSQW = 0;
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }

    private static final long serialVersionUID = 1L;
    protected Space space;
    protected IntegratorHard integratorHard;
    protected double virialSumIG, virialSumSQW;
    protected double lastTime;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
}
