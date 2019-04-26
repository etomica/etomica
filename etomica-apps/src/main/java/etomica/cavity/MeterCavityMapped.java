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
    protected final Vector deltaMomentum;
    protected final IntegratorHard integratorHard;
    protected double lastTime;
    protected final Vector dr;
    protected double[] gSum;
    protected double lastSwitchTime;
    protected boolean internal;
    protected double tInternal, tExternal;
    protected double sigma;
    protected double mappingCut2 = Double.POSITIVE_INFINITY;
    protected boolean resetAfterData;
    protected long internalCollisions, totalCollisions;
    protected IDataSink rawSink;
    protected final DataTag rawTag;
    protected DataFunction.DataInfoFunction rawDataInfo;
    protected DataDoubleArray rawData;
    protected boolean useMomentum = false;

    public MeterCavityMapped(IntegratorHard integrator) {

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setNValues(500);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.INCLUSIVE);
        xDataSource.setTypeMin(DataSourceUniform.LimitType.INCLUSIVE);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new double[rData.getLength()];

        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("mapped y(r)", Null.DIMENSION, this);
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
        dataInfo = new DataFunction.DataInfoFunction("mapped y(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        if (rawSink != null) {
            setRawSink(rawSink);
        }
        zeroData();
    }

    public void zeroData() {
        for (int i = 0; i < gSum.length; i++) gSum[i] = 0;
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
        P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;
        P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
        sigma = p2.getCollisionDiameter();
        if (cType == P2HardSphereCavity.CollisionType.CAPTURE) {
            internal = true;
            double t = integratorHard.getCurrentTime() + agent.collisionTime();
            tExternal += t - lastSwitchTime;
            lastSwitchTime = t;
        }
        else if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
            internalCollisions++;
            internal = false;
            double t = integratorHard.getCurrentTime() + agent.collisionTime();
            tInternal += t - lastSwitchTime;
            lastSwitchTime = t;
        } else if (cType == P2HardSphereCavity.CollisionType.INTERNAL_BOUNCE) {
            internalCollisions++;
        }

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

        if (useMomentum) {
            deltaMomentum.Ea1Tv1(p2.lastCollisionVirial() / dr.squared(), dr);
        } else {
            deltaMomentum.Ea1Tv1(-1.0 / sigma, dr);
        }

        if (atom1Paired) {
            if (atom1 == p2.pairedAtom1) {
                atom2 = p2.pairedAtom2;
            } else {
                atom2 = p2.pairedAtom1;
            }
        } else {
            // fix virial to be for atom2
            deltaMomentum.TE(-1);
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
        if (r2 > mappingCut2) return;
        double r = Math.sqrt(r2);

        int index = (int) (r / sigma * (xDataSource.getNValues() - 1));
        if (atom2Paired) dr.TE(-1);
        gSum[index] += deltaMomentum.dot(dr) / (r * r2);
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
        int N = integratorHard.getBox().getLeafList().getAtoms().size();
        double V = integratorHard.getBox().getBoundary().volume();
        IData rData = xDataSource.getData();
        double dx = rData.getValue(2) - rData.getValue(1);
        double density = N / V;
        double sqrtPI = Math.sqrt(Math.PI);
        for (int i = 0; i < y.length; i++) {
            if (useMomentum) {
                // also divide by T
                y[i] = gSum[i] / elapsedTime / (N * density * 4 * Math.PI);
            } else {
                // also multiply by (T/mass)^.5
                y[i] = gSum[i] / elapsedTime / (N * density * 4 * sqrtPI);
            }
        }
        long externalCollision = totalCollisions - internalCollisions;
        double fac = externalCollision / (double) internalCollisions;
        if (rawSink != null) {
            double[] f = rawData.getData();
            for (int i = 0; i < f.length; i++) {
                double r = rData.getValue(i) + dx * 0.5;
                f[i] = y[i];
                if (r < sigma) f[i] *= fac;
            }
            rawSink.putData(rawData);
        }

        for (int i = 0; i < y.length; i++) {
            for (int j = i + 1; j < y.length; j++) {
                y[i] += y[j];
            }
        }
        double yIntegral = 0;
        for (int i = 0; i < y.length; i++) {
            double r = rData.getValue(i);
            yIntegral += r * r * y[i];
        }

        yIntegral *= 4 * Math.PI * dx;
        double nPairs = yIntegral * N * (N - 1) / V / 2;
//        System.out.println("yIntegral "+yIntegral+" nPairs: "+nPairs);
        // need to shift
        double ti = tInternal, te = tExternal;
        if (internal) ti += integratorHard.getCurrentTime() - lastSwitchTime;
        else te += integratorHard.getCurrentTime() - lastSwitchTime;
        double shift = ((ti / (te + ti)) - nPairs) / (4.0 / 3.0 * Math.PI * N * (N - 1) / V * sigma * sigma * sigma / 2);
        if (ti * te > 0) {
//            System.out.println("shift by "+shift+" so nPairs "+nPairs+" matches "+(ti/(te+ti))+" "+ti+" "+te);
            for (int i = 0; i < y.length; i++) {
                y[i] += shift;
            }
        }
        if (false) {
            // recheck integral
            yIntegral = 0;
            // the last point is our endpoint.  skip and handle explicitly
            for (int i = 0; i < y.length - 1; i++) {
                double r = rData.getValue(i);
                yIntegral += r * r * y[i];
            }
            // our end point isn't 0 anymore
            yIntegral += 0.5 * shift * sigma * sigma;

            yIntegral *= 4 * Math.PI * dx;
            nPairs = yIntegral * N * density / 2;
            System.out.println("and we got integral " + yIntegral + " and nPairs " + nPairs);
        }
        // and now scale
        if (internalCollisions > 0) {

            for (int i = 0; i < y.length; i++) {
                y[i] *= fac;
            }
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
