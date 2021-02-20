/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorHard;
import etomica.potential.P2HardSphere;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for the cavity function using mapped averaging.
 *
 * @author Andrew Schultz
 */
public class MeterRDFMapped implements IDataSource, IntegratorHard.CollisionListener, DataSourceIndependent {

    protected final DataSourceUniform xDataSource;
    protected DataDoubleArray rData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataFunction data;
    protected DataTag tag;
    protected final Vector deltaMomentum;
    protected final IntegratorHard integratorHard;
    protected double lastTime;
    protected final Vector dr;
    protected double[] gSum, gSum2;
    protected double lastSwitchTime;
    protected boolean internal;
    protected double tInternal, tExternal;
    protected double sigma;
    protected boolean resetAfterData;
    protected long internalCollisions, totalCollisions;
    protected IDataSink rawSink;
    protected final DataTag rawTag;
    protected DataFunction.DataInfoFunction rawDataInfo;
    protected DataDoubleArray rawData;
    public boolean foobar = false;
    protected double mappingCut2 = Double.POSITIVE_INFINITY;
    protected boolean useMomentum;

    public MeterRDFMapped(IntegratorHard integrator) {

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setNValues(500);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.INCLUSIVE);
        xDataSource.setTypeMin(DataSourceUniform.LimitType.INCLUSIVE);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new double[rData.getLength()];
        gSum2 = new double[rData.getLength()];

        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("mapped g(r)", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        lastTime = integratorHard.getCurrentTime();
        dr = integrator.getBox().getSpace().makeVector();
        deltaMomentum = integrator.getBox().getSpace().makeVector();
        rawTag = new DataTag();
    }

    public void setRawSink(IDataSink rawSink) {
        this.rawSink = rawSink;
        if (rawSink == null) {
            rawData = null;
            rawDataInfo = null;
            return;
        }
        rawData = new DataFunction(new int[]{rData.getLength()});
        rawDataInfo = new DataFunction.DataInfoFunction("f(r)", Null.DIMENSION, this);
        rawDataInfo.addTag(rawTag);
        rawSink.putDataInfo(rawDataInfo);
    }

    /**
     * Resets the data collected.
     */
    public void reset() {

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new double[rData.getLength()];
        gSum2 = new double[rData.getLength()];
        dataInfo = new DataFunction.DataInfoFunction("mapped g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        if (rawSink != null) {
            setRawSink(rawSink);
        }
        zeroData();
    }

    public void zeroData() {
        for (int i = 0; i < gSum.length; i++) gSum[i] = gSum2[i] = 0;
        lastTime = integratorHard.getCurrentTime();
        lastSwitchTime = integratorHard.getCurrentTime();
        tInternal = 0;
        tExternal = 0;
        internalCollisions = totalCollisions = 0;
    }

    public void setMappingCut(double mappingCut) {
        mappingCut2 = mappingCut * mappingCut;
    }

    public void setResetAfterData(boolean doResetAfterData) {
        resetAfterData = doResetAfterData;
    }

    public void collisionAction(IntegratorHard.Agent agent) {
        totalCollisions++;
        IAtomKinetic atom1 = agent.atom;
        IAtomKinetic atom2 = agent.collisionPartner;
        IAtomKinetic[] atoms12 = new IAtomKinetic[]{atom1, atom2};
        sigma = ((P2HardSphere) agent.collisionPotential).getCollisionDiameter();

        if (agent.collisionPotential instanceof P2HardSphereCavity) {
            P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;
            P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
            sigma = p2.getCollisionDiameter();
            if (cType == P2HardSphereCavity.CollisionType.CAPTURE) {
                internal = true;
                double t = integratorHard.getCurrentTime() + agent.collisionTime();
                tExternal += t - lastSwitchTime;
                lastSwitchTime = t;
            } else if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
                internalCollisions++;
                internal = false;
                double t = integratorHard.getCurrentTime() + agent.collisionTime();
                tInternal += t - lastSwitchTime;
                lastSwitchTime = t;
            } else if (cType == P2HardSphereCavity.CollisionType.INTERNAL_BOUNCE) {
                internalCollisions++;
            }
        }

        double falseTime = agent.collisionTime();

        Box box = integratorHard.getBox();
        Boundary boundary = box.getBoundary();

        // actual colliding pair.  virial is for atom1
        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        dr.PEa1Tv1(falseTime, atom1.getVelocity());
        dr.PEa1Tv1(-falseTime, atom2.getVelocity());
        boundary.nearestImage(dr);

        if (useMomentum) {
            deltaMomentum.Ea1Tv1(agent.collisionPotential.lastCollisionVirial() / dr.squared(), dr);
        } else {
            deltaMomentum.Ea1Tv1(-1.0 / sigma, dr);
        }

        Vector r12 = null;
        if (foobar) {
            r12 = new Vector3D();
            r12.E(atom2.getPosition());
            r12.PEa1Tv1(falseTime, atom2.getVelocity());
            r12.PEa1Tv1(0.5, dr);
        }
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms12.length; i++) {
            for (int j = 0; j < atoms.size(); j++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(j);
                if (a == atom1 || a == atom2) continue;
                dr.Ev1Mv2(atoms12[i].getPosition(), a.getPosition());
                dr.PEa1Tv1(falseTime, atoms12[i].getVelocity());
                dr.PEa1Tv1(-falseTime, a.getVelocity());
                boundary.nearestImage(dr);

                double r2 = dr.squared();
                if (r2 > mappingCut2) continue;
                double r = Math.sqrt(r2);
                int index = (int) (r / xDataSource.getXMax() * (xDataSource.getNValues() - 1));
                if (index >= gSum.length) index = gSum.length - 1;
                double x = deltaMomentum.dot(dr) / (r * r2);
                gSum[index] += x;

                if (foobar) {
                    dr.Ev1Mv2(a.getPosition(), r12);
                    dr.PEa1Tv1(falseTime, a.getVelocity());
                    boundary.nearestImage(dr);
                    double foor2 = dr.squared();
                    double foor = Math.sqrt(foor2);
                    index = (int) (foor / xDataSource.getXMax() * (xDataSource.getNValues() - 1));
                    if (index >= gSum.length) index = gSum.length - 1;
                    gSum2[index] += x;
                }
            }
            deltaMomentum.TE(-1);
        }
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
        int N = integratorHard.getBox().getLeafList().size();
        double V = integratorHard.getBox().getBoundary().volume();
        IData rData = xDataSource.getData();
        double dx = rData.getValue(2) - rData.getValue(1);
        double density = N / V;
        double[] f = rawData == null ? null : rawData.getData();
        double sqrtPI = Math.sqrt(Math.PI);
        for (int i = 0; i < y.length; i++) {
            y[i] = gSum[i] / (elapsedTime * N * density * 4);
            if (foobar) f[i] = gSum2[i] / (elapsedTime * N * density * 4);
            if (useMomentum) {
                // also divide by T
                y[i] /= Math.PI;
                if (foobar) f[i] /= Math.PI;
            } else {
                // also multiply by (T/mass)^.5
                y[i] /= sqrtPI;
                if (foobar) f[i] /= sqrtPI;
            }
        }

        if (rawSink != null) {
            if (!foobar) {
                for (int i = 0; i < f.length; i++) {
                    f[i] = y[i];
                }
            }
            rawSink.putData(rawData);
        }

        long externalCollision = totalCollisions - internalCollisions;
        double fac = externalCollision / (double) internalCollisions;

        if (internalCollisions > 0) {
            for (int i = 0; i < y.length; i++) {
                double r = rData.getValue(i) + dx * 0.5;
                if (r < sigma) y[i] *= fac;
            }
        }

        for (int i = 0; i < y.length; i++) {
            for (int j = i + 1; j < y.length; j++) {
                y[i] += y[j];
            }
            y[i] += 1;
        }
        if (resetAfterData) zeroData();
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
