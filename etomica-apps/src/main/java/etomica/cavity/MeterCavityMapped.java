/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for the cavity function using mapped averaging.
 *
 * @author Andrew Schultz
 */
public class MeterCavityMapped implements IDataSource, IntegratorHard.CollisionListener, DataSourceIndependent {

    protected final DataSourceUniform xDataSource;
    protected DataDoubleArray rData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataFunction data;
    protected DataTag tag;
    protected final Vector lastCollisionVirial;
    protected final IntegratorHard integratorHard;
    protected double lastTime;
    protected final Vector dr;
    protected double[] gSum;

    public MeterCavityMapped(IntegratorHard integrator) {

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.HALF_STEP);
        xDataSource.setTypeMin(DataSourceUniform.LimitType.HALF_STEP);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new double[rData.getLength()];

        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("cavity function", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        lastTime = integratorHard.getCurrentTime();
        dr = integrator.getBox().getSpace().makeVector();
        lastCollisionVirial = integrator.getBox().getSpace().makeVector();
    }

    /**
     * Resets the data collected.
     */
    public void reset() {

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new double[rData.getLength()];
        dataInfo = new DataFunction.DataInfoFunction("mapped cavity(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);

        lastTime = integratorHard.getCurrentTime();
    }

    public void collisionAction(IntegratorHard.Agent agent) {
        IAtomKinetic atom1 = agent.atom;
        IAtomKinetic atom2 = agent.collisionPartner;
        P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;

        boolean atom1Paired = atom1 == p2.pairedAtom1 || atom1 == p2.pairedAtom2;
        boolean atom2Paired = atom2 == p2.pairedAtom1 || atom2 == p2.pairedAtom2;
        if (atom1Paired == atom2Paired) return;
        double falseTime = agent.collisionTime();

        // actual colliding pair.  virial is for atom1
        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        dr.PEa1Tv1(falseTime, atom1.getVelocity());
        dr.PEa1Tv1(-falseTime, atom2.getVelocity());
        Box box = integratorHard.getBox();
        box.getBoundary().nearestImage(dr);
        lastCollisionVirial.Ea1Tv1(p2.lastCollisionVirial() / dr.squared(), dr);

        if (atom1Paired) {
            if (atom1 == p2.pairedAtom1) {
                atom2 = p2.pairedAtom2;
            } else {
                atom2 = p2.pairedAtom1;
            }
        } else {
            // fix virial to be for atom2
            lastCollisionVirial.TE(-1);
            if (atom2 == p2.pairedAtom1) {
                atom1 = p2.pairedAtom2;
            } else {
                atom1 = p2.pairedAtom1;
            }
        }
        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        dr.PEa1Tv1(falseTime, atom1.getVelocity());
        dr.PEa1Tv1(-falseTime, atom2.getVelocity());
        integratorHard.getBox().getBoundary().nearestImage(dr);
        double r2 = dr.squared();
        double r = Math.sqrt(r2);
        int index = xDataSource.getIndex(r);  //determine histogram index
        if (atom2Paired) dr.TE(-1);
        gSum[index] += lastCollisionVirial.dot(dr) / (r * r2);
    }

    public IData getData() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength()) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }

        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastTime;
        if (elapsedTime == 0.0) return data;
        if (elapsedTime < 0) throw new RuntimeException("you should have called reset");

        final double[] y = data.getData();
        double[] r = rData.getData();
        double dx2 = 0.5 * (xDataSource.getXMax() - xDataSource.getXMin()) / r.length;
        Space space = integratorHard.getBox().getSpace();
        int N = integratorHard.getBox().getLeafList().getAtoms().size();
        for (int i = 0; i < r.length; i++) {
            double vShell = space.sphereVolume(r[i] + dx2) - space.sphereVolume(r[i] - dx2);
            y[i] = gSum[i] / (vShell * elapsedTime) / N;
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

}
