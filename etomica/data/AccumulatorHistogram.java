/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica.data;

import etomica.DataTranslator;
import etomica.utility.Histogram;
import etomica.utility.HistogramSimple;

/**
 * Accumulator that keeps histogram of data.
 */
public class AccumulatorHistogram extends DataAccumulator {
	
	Histogram[] histogram = new Histogram[0];
	private double[][] data;
	int nData, nDataMinus1;
	private String xLabel;
	private Histogram.Factory histogramFactory;
	private int nBins;
	private DataTranslatorArray dataTranslator;
	
	/**
	 * Creates instance using HistogramSimple factory and specifying
	 * histograms having 100 bins.
	 */
	public AccumulatorHistogram() {
		this(HistogramSimple.FACTORY);
	}
	public AccumulatorHistogram(Histogram.Factory factory) {
		this(factory, 100);
	}
	public AccumulatorHistogram(Histogram.Factory factory, int nBins) {
		this.nBins = nBins;
		setNData(0);
		histogramFactory = factory;
	}

	/**
	 * Adds each value in the given array to its own histogram.
	 * If the number of values is different from that given in 
	 * previous calls to the method, old histogram data is discarded
	 * and new histograms are constructed (this behavior can be modified
	 * by overriding the setNData method).
	 */
	protected void addData(double[] values) {
		if(values.length != nData) setNData(values.length);
		for(int i=nDataMinus1; i>=0; i--) {       		
			histogram[i].addValue(values[i]);
		}
	}

	/**
	 * Returns the set of histograms consecutively in a 1D array.
	 */
	public double[] getData() {
		for(int i=0; i<nData; i++) data[i] = histogram[i].getHistogram();
		return dataTranslator.toArray(data);
	}

	public void setLabel(String s) {label = s;}
	public String getLabel() {return label;}
	       
	public void setXLabel(String s) {xLabel = s;}
	public String getXLabel() {return xLabel;}
        

	/**
	 * Returns a DataTranslatorArray instance with dimensions
	 * appropriate to current values of nData and nBins.
	 */
	public DataTranslator getTranslator() {
		return new DataTranslatorArray(nData, nBins);
	}
	
	/**
	 * Determines the number of histograms to be recorded, and is
	 * invoked when the add method is called with an array argument
	 * of length different than previously given.  As implemented here,
	 * setNData(int) causes all previously-recorded histograms to
	 * be discarded.
	 * @param nData
	 */
	protected void setNData(int nData) {
    	this.nData = nData;
    	nDataMinus1 = nData-1;
    	dataTranslator = (DataTranslatorArray)getTranslator();
    	data = new double[nData][];
    	histogram = new Histogram[nData];
    	for(int i=0; i<nData; i++) histogram[i] = histogramFactory.makeHistogram(nBins);
	}
	/**
	 * @return Returns nBins, the number of bins in each histogram.
	 */
	public int getNBins() {
		return nBins;
	}
	/**
	 * @param bins Sets the number of bins in each histogram.  Calls
	 * setNBins method of current histograms, which will discard data or
	 * modify themselves depending on how they are defined. 
	 */
	public void setNBins(int nBins) {
		this.nBins = nBins;
    	dataTranslator = (DataTranslatorArray)getTranslator();
		for(int i=0; i<nData; i++) histogram[i].setNBins(nBins);
	}
	
	public void reset() {
		for (int i=0; i<nData; i++) histogram[i].reset();
	}
}
