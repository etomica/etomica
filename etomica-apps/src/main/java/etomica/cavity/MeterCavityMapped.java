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
    protected P2HardSphereCavity p2;

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
        deltaMomentum = integrator.getBox().getSpace().makeVector();
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
        tInternal = 0;
        tExternal = 0;

    }

    public void collisionAction(IntegratorHard.Agent agent) {
        IAtomKinetic atom1 = agent.atom;
        IAtomKinetic atom2 = agent.collisionPartner;
        p2 = (P2HardSphereCavity) agent.collisionPotential;
        P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
        sigma = p2.getCollisionDiameter();
        if (cType == P2HardSphereCavity.CollisionType.CAPTURE) {
            internal = true;
            tExternal += integratorHard.getCurrentTime() - lastSwitchTime;
        }
        else if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
            internal = false;
            tInternal += integratorHard.getCurrentTime() - lastSwitchTime;
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

        deltaMomentum.Ea1Tv1(p2.lastCollisionVirial() / dr.squared(), dr);

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
        double r = Math.sqrt(r2);
        int index = xDataSource.getIndex(r);  //determine histogram index
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
        for (int i = 0; i < y.length; i++) {
            y[i] = gSum[i] / elapsedTime / (N * density * 4 * Math.PI);
        }

        for (int i = 0; i < y.length; i++) {
            for (int j = i + 1; j < y.length; j++) {
                y[i] += y[j];
            }
        }
//        if (true) return data;
        double yIntegral = 0;
        for (int i = 0; i < y.length; i++) {
            double r = rData.getValue(i);
            yIntegral += r * r * y[i];
        }

        yIntegral *= 4 * Math.PI * N * density / 3 * dx / 2;
//        System.out.println("yIntegral "+yIntegral);
        // need to shift
        double ti = tInternal, te = tExternal;
        if (internal) ti += integratorHard.getCurrentTime() - lastSwitchTime;
        else te += integratorHard.getCurrentTime() - lastSwitchTime;
        if (ti * te > 0) {
//            System.out.println("shifting so integral is "+ti/(te+ti));
            double shift = ((ti / (te + ti)) - yIntegral) / (4 * Math.PI * N * density * sigma * sigma * sigma / 3 / 2);
//            System.out.println("gonna shift by "+shift);
            for (int i = 0; i < y.length; i++) {
                y[i] += shift;
            }
        }
//        yIntegral = 0;
//        for (int i = 0; i < y.length; i++) {
//            double r = rData.getValue(i);
//            yIntegral += r*r*y[i];
//        }
//        yIntegral *= 4*Math.PI*N*density/3*dx/2;
//        System.out.println("yIntegral => "+yIntegral);
        if (p2 != null) {
            long totalCollision = integratorHard.getCollisionCount();
            long internalCollision = p2.getInternalCount();
            double fac = 1;
            if (internalCollision > 0) {
                long externalCollision = totalCollision - internalCollision;
                fac = externalCollision / (double) internalCollision;
            }
            for (int i = 0; i < y.length; i++) {
                y[i] *= fac;
            }
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
