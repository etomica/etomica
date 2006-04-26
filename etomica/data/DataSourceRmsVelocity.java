package etomica.data;

import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.ICoordinateKinetic;
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

public class DataSourceRmsVelocity implements DataSourceAtomic, Serializable {

    public DataSourceRmsVelocity() {
        this(new HistogramCollapsing());
    }
    
	public DataSourceRmsVelocity(Histogram histogram) {
        setIterator(null);
        atomData = new DataDouble("RMS Velocity", new DimensionRatio(Length.DIMENSION, Time.DIMENSION));
        this.histogramRMS = histogram;
        setupData();
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

	/**
	 * Returns the rms velocity of the atoms given by the iterator. Value is
	 * given in the first element of the array, which is always of dimension 1.
	 */
	public Data getData() {
		iterator.reset();
        histogramRMS.reset();
		while (iterator.hasNext()) {
			AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
			histogramRMS.addValue(Math.sqrt(((ICoordinateKinetic)atom.coord).velocity().squared()));
		}

        //covertly invoke getHistogram, which actually calculates the histogram
        // our DataFunction wraps the histogram array
        if (data.getData() != histogramRMS.getHistogram() ||
                data.getXData(0).getData() != histogramRMS.xValues()) {
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
        DataDoubleArray xData = new DataDoubleArray(atomData.getDataInfo().getLabel(),atomData.getDataInfo().getDimension(), 
                                                    new int[]{nBins},histogramRMS.xValues());
        data = new DataFunction("RMS Velocity Histogram",Null.DIMENSION, new DataDoubleArray[]{xData}, histogramRMS.getHistogram());
    }
    
    public Data getData(Atom a) {
        atomData.x = Math.sqrt(((ICoordinateKinetic)((AtomLeaf)a).coord).velocity().squared());
        return atomData;
    }        
    
	/**
	 * Sets the iterator defining the atoms for which the RMS velocity is
	 * calculated.
	 * 
	 * @param iterator
	 */
	public void setIterator(AtomIterator iterator) {
		if (iterator == null)
			this.iterator = AtomIteratorNull.INSTANCE;
		else
			this.iterator = iterator;
	}

	/**
	 * @return The iterator defining the atoms used for the RMS velocity
	 *         calculation
	 */
	public AtomIterator getIterator() {
		return iterator;
	}

	private AtomIterator iterator;
    private final DataDouble atomData;
    private final Histogram histogramRMS;
    private DataFunction data;

}//end of DataSourceVelocityRms
