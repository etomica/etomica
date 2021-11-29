/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  This
 * implementation uses a NeighborManager that must be configured to loop over
 * neighbors within the range of the RDF.
 */
public class MeterRDFNeighbors implements IDataSource, DataSourceIndependent {

    /**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a box.
	 */
    public MeterRDFNeighbors(Box box, NeighborManager neighborManager) {
        this.box = box;
	    this.neighborManager = neighborManager;
        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.HALF_STEP);
        xDataSource.setTypeMin(DataSourceUniform.LimitType.HALF_STEP);

        gSum = new long[xDataSource.getData().getLength()];

        rData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {rData.getLength()});
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

        tag = new DataTag();
        dataInfo.addTag(tag);
        neighborIterator = neighborManager.makeNeighborIterator();
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
        gSum = new long[xDataSource.getData().getLength()];
        rData = (DataDoubleArray)xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[] {rData.getLength()});
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        callCount = 0;
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
        }

        int numAtoms = box.getLeafList().size();
        double xMaxSquared = xMax*xMax;
        for (int i = 0; i<numAtoms; i++) {
            neighborIterator.iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    double r2 = rij.squared();       //compute pair separation
                    if(r2 < xMaxSquared) {
                        int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                        gSum[index]++;                        //add once for each atom
                    }
                }
            });
        }

        callCount++;

        final double[] y = data.getData();
        long numAtomPairs = numAtoms*(numAtoms-1)/2;
	    double norm = numAtomPairs * callCount / box.getBoundary().volume();
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    Space space = box.getSpace();
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

    protected long[] gSum;
    protected Box box;
    protected DataFunction data;
    private IDataInfo dataInfo;
    protected DataDoubleArray rData;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    protected final DataTag tag;
    protected long callCount;
    protected final NeighborManager neighborManager;
    protected final NeighborIterator neighborIterator;
}
