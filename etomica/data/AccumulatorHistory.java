/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica.data;

import etomica.DataSource;
import etomica.DataTranslator;
import etomica.utility.History;
import etomica.utility.HistoryScrolling;

/**
 * Accumulator that keeps history of data.
 */
public class AccumulatorHistory extends DataAccumulator {
	
	History[] history = new History[0];
    DataSourceUniform xSource;
	private double[][] data;
	int nData, nDataMinus1;
	private String xLabel;
	private History.Factory historyFactory;
	private int historyLength;
	private DataTranslatorArray dataTranslator;
	
	/**
	 * Creates instance using HistorySimple factory and specifying
	 * historys having 100 bins.
	 */
	public AccumulatorHistory() {
		this(HistoryScrolling.FACTORY);
	}
	public AccumulatorHistory(History.Factory factory) {
		this(factory, 100);
	}
	public AccumulatorHistory(History.Factory factory, int historyLength) {
	    super();
        xSource = new DataSourceUniform(1,historyLength,historyLength);
		this.historyLength = historyLength;
		historyFactory = factory;
        setNData(1);
	}

	/**
	 * Adds each value in the given array to its own history.
	 * If the number of values is different from that given in 
	 * previous calls to the method, old history data is discarded
	 * and new historys are constructed (this behavior can be modified
	 * by overriding the setNData method).
	 */
	protected void addData(double[] values) {
		if(values.length != nData) setNData(values.length);
		for(int i=nDataMinus1; i>=0; i--) {       		
			history[i].addValue(values[i]);
		}
	}

	/**
	 * Returns the set of histories consecutively in a 1D array.
	 */
	public double[] getData() {
		for(int i=0; i<nData; i++) data[i] = history[i].getHistory();
		return dataTranslator.toArray(data);
	}

	public void setLabel(String s) {label = s;}
	public String getLabel() {return label;}
	       
	public void setXLabel(String s) {xLabel = s;}
	public String getXLabel() {return xLabel;}
    
    public DataSource getXSource() {
        return history[0].getXSource();
    }

	/**
	 * Returns a DataTranslatorArray instance with dimensions
	 * appropriate to current values of nData and historyLength.
	 */
	public DataTranslator getTranslator() {
		return new DataTranslatorArray(nData, historyLength);
	}
	
	/**
	 * Determines the number of historys to be recorded, and is
	 * invoked when the add method is called with an array argument
	 * of length different than previously given.  As implemented here,
	 * setNData(int) causes all previously-recorded histories to
	 * be discarded.
	 * @param nData
	 */
	protected void setNData(int nData) {
    	this.nData = nData;
    	nDataMinus1 = nData-1;
    	dataTranslator = (DataTranslatorArray)getTranslator();
    	data = new double[nData][];
    	history = new History[nData];
    	for(int i=0; i<nData; i++) history[i] = historyFactory.makeHistory(historyLength);
	}
	/**
	 * @return Returns historyLength, the number of bins in each history.
	 */
	public int getHistoryLength() {
		return historyLength;
	}
	/**
	 * @param bins Sets the number of bins in each history.  Calls
	 * setHistoryLength method of current histories, which will discard data or
	 * modify themselves depending on how they are defined. 
	 */
	public void setHistoryLength(int historyLength) {
		this.historyLength = historyLength;
        xSource.setNValues(historyLength);
        xSource.setXMax(historyLength);
    	dataTranslator = (DataTranslatorArray)getTranslator();
		for(int i=0; i<nData; i++) history[i].setHistoryLength(historyLength);
	}
	
	public void reset() {
    	for(int i=0; i<nData; i++) history[i].reset();
	}
}
