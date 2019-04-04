/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorBoxDependent;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for tabulating the cavity function for hard spheres based on tracking
 * how often various pair separations are visited.
 *
 * @author Andrew Schultz
 */
public class MeterCavity implements IAction, IDataSource, DataSourceIndependent {

    protected final Space space;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr;
    protected final Box box;
    protected long[] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected AtomsetIteratorBoxDependent iterator;
    protected double xMax;
    protected long callCount;
    private IDataInfo dataInfo;
    private final P2HardSphereCavity p2;

    /**
     * Creates meter with default to compute pair correlation for all
     * leaf atoms in a box.
     */
    public MeterCavity(Box box, P2HardSphereCavity p2) {
        this.box = box;
        this.space = box.getSpace();
        this.p2 = p2;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        xDataSource.setXMax(p2.getCollisionDiameter());

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

        iterator = new ApiLeafAtoms();
        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray) xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("cavity(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        callCount = 0;
    }

    /**
     * Takes the RDF for the current configuration of the given box.
     */
    public void actionPerformed() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
        }

        IAtom atom1 = p2.getPairedAtom1();
        if (atom1 == null) return;
        IAtom atom2 = p2.getPairedAtom2();
        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        box.getBoundary().nearestImage(dr);
        double r2 = dr.squared();
        double xMax = p2.getCollisionDiameter();
        if (r2 < xMax * xMax) {
            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
            gSum[index]++;                        //add once for each atom
        }
        callCount++;
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }

        final double[] y = data.getData();
        long numAtomPairs = 1;
        double sigma = p2.getCollisionDiameter();
        double vol = 4.0 / 3.0 * sigma * sigma * sigma;
        double norm = numAtomPairs * callCount / vol;
        double[] r = rData.getData();
        double dx2 = 0.5 * (xMax - xDataSource.getXMin()) / r.length;
        for (int i = 0; i < r.length; i++) {
            double vShell = space.sphereVolume(r[i] + dx2) - space.sphereVolume(r[i] - dx2);
            y[i] = gSum[i] / (norm * vShell);
        }
        return data;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
}
