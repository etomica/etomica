/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.action.IAction;
import etomica.api.IBoundary;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
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
import etomica.units.Length;
import etomica.units.Null;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.
 *
 * @author David Kofke
 */
public class MeterRDF implements IAction, IEtomicaDataSource, DataSourceIndependent, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    private final Vector dr;
    protected Box box;
    protected long[] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected AtomsetIteratorBoxDependent iterator;
    protected double xMax;
    protected long callCount;
    protected AtomType type1, type2;
    private IEtomicaDataInfo dataInfo;
    private IBoundary boundary;
    private String name;

	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a box.
	 * @param space
	 */
    public MeterRDF(Space space) {
	    this.space = space;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        rData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

	    iterator = new ApiLeafAtoms();
        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setAtomType(AtomType type) {
        type1 = type;
        type2 = type;
    }

    public void setAtomTypes(AtomType type1, AtomType type2) {
        this.type1 = type1;
        this.type2 = type2;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray)xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[] {rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
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

        double xMaxSquared = xMax*xMax;
        iterator.setBox(box);
        iterator.reset();
        // iterate over all pairs
        for (IAtomList pair = iterator.next(); pair != null;
             pair = iterator.next()) {
            if (type1 != null && (pair.getAtom(0).getType() != type1 || pair.getAtom(1).getType() != type2)) continue;
            dr.Ev1Mv2(pair.getAtom(1).getPosition(),pair.getAtom(0).getPosition());
            boundary.nearestImage(dr);
            double r2 = dr.squared();       //compute pair separation
            if(r2 < xMaxSquared) {
                int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                gSum[index]++;                        //add once for each atom
            }
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
        long numAtomPairs = 0;
        if (type1 == null) {
            long numAtoms = box.getLeafList().getAtomCount();
            numAtomPairs = numAtoms*(numAtoms-1)/2;
        }
        else {
            iterator.setBox(box);
            iterator.reset();
            for (IAtomList pair = iterator.next(); pair != null; pair = iterator.next()) {
                if (pair.getAtom(0).getType() != type1 || pair.getAtom(1).getType() != type2) continue;
                numAtomPairs++;
            }
        }
	    double norm = numAtomPairs * callCount / box.getBoundary().volume();
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    for(int i=0;i<r.length; i++) {
	        double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
	        y[i] = gSum[i] / (norm*vShell);
	    }
	    return data;
	}

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
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

    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
