/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IAtomKinetic;
import etomica.atom.iterator.AtomIterator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.DimensionRatio;
import etomica.units.Length;
import etomica.units.Null;
import etomica.units.Time;
import etomica.util.Histogram;
import etomica.util.HistogramCollapsing;

/**
 * Meter for the root-mean-square velocity of a set of atoms. Useful to obtain
 * histograms of atom speeds.
 * 
 * @author David Kofke
 */

public class DataSourceRmsVelocity implements IEtomicaDataSource, DataSourceAtomic, DataSourceIndependent, Serializable {

    public DataSourceRmsVelocity() {
        this(new HistogramCollapsing());
    }
    
	public DataSourceRmsVelocity(Histogram histogram) {
        setIterator(null);
        atomDataInfo = new DataInfoDouble("RMS Velocity", new DimensionRatio(Length.DIMENSION, Time.DIMENSION));
        atomData = new DataDouble();
        this.histogramRMS = histogram;
        tag = new DataTag();
        atomDataInfo.addTag(tag);
        xTag = new DataTag();
        setupData();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getAtomDataInfo() {
        return atomDataInfo;
    }

	/**
	 * Returns the rms velocity of the atoms given by the iterator. Value is
	 * given in the first element of the array, which is always of dimension 1.
	 */
	public IData getData() {
		iterator.reset();
        histogramRMS.reset();
		for (IAtomKinetic atom = (IAtomKinetic)iterator.nextAtom(); atom != null;
             atom = (IAtomKinetic)iterator.nextAtom()) {
			histogramRMS.addValue(Math.sqrt(atom.getVelocity().squared()));
		}

        //covertly invoke getHistogram, which actually calculates the histogram
        // our DataFunction wraps the histogram array
        if (data.getData() != histogramRMS.getHistogram() ||
                xData.getData() != histogramRMS.xValues()) {
            // we wrap the histogram inner array instances in the DataFunction.
            // if they change, we need a new instance of the DataFunction
            setupData();
        }
        
		return data;
	}
    
    /**
     * Creates the data object (a DataFunction) to be returned by getData().  
     * data wraps the histogram's double[] so copying is not needed.
     */
    protected void setupData() {
        int nBins = histogramRMS.getNBins();
        xData = new DataDoubleArray(new int[]{nBins},histogramRMS.xValues());
        xDataInfo = new DataInfoDoubleArray(atomDataInfo.getLabel(),atomDataInfo.getDimension(), new int[]{nBins});
        xDataInfo.addTag(xTag);
        data = new DataFunction(new int[]{nBins}, histogramRMS.getHistogram());
        dataInfo = new DataInfoFunction("RMS Velocity Histogram",Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }
    
    public IData getData(IAtom a) {
        atomData.x = Math.sqrt(((IAtomKinetic)a).getVelocity().squared());
        return atomData;
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

	/**
	 * Sets the iterator defining the atoms for which the RMS velocity is
	 * calculated.
	 * 
	 * @param newIterator
	 */
	public void setIterator(AtomIterator newIterator) {
	    iterator = newIterator;
	}

	/**
	 * @return The iterator defining the atoms used for the RMS velocity
	 *         calculation
	 */
	public AtomIterator getIterator() {
		return iterator;
	}

    private static final long serialVersionUID = 1L;
	private AtomIterator iterator;
    private final DataDouble atomData;
    private final IEtomicaDataInfo atomDataInfo;
    private IEtomicaDataInfo dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    private final Histogram histogramRMS;
    private DataFunction data;
    protected final DataTag tag, xTag;
}
