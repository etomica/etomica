/*
 * History
 * Created on Aug 9, 2004 by kofke
 */
package etomica.utility;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface History {
	public void setNValues(int n);

	public int getNValues();

	public void reset();

	public void addValue(double x);

	public double[] xValues();

	public double[] getHistory();
	
	public interface Factory {
		public History makeHistory();
		public History makeHistory(int n);
	}

}