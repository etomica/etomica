/*
 * Created on Jul 22, 2003
 */
package etomica;

import etomica.units.Dimension;
import etomica.utility.*;
/**
 *
 * Secondary Meter which monitors changes of the first transform over time.
 * It then transforms this data to result in the frequency
 * 
 * @author Mitchell Tai
 */
public class MeterFourierTransform2 extends MeterGroup {
	
	private MeterFourierTransform transform;
	private int nMeters;
	private double[] current;
	
	/**
	 * @param sim
	 * @param nMeters
	 * @param transform
	 */
	public MeterFourierTransform2(Simulation sim, int nMeters,MeterFourierTransform transform) {
		super(sim, nMeters);
		setResetHistoryOnMeterReset(true);
		setHistorying(true);
		for (int i=0;i<nMeters;i++) {
			getHistory(i).setNValues(2*nMeters);
			getHistory(i).setLabel("History " + (i+1));
			getHistory(i).setTransform(new Transform.Fourier());
			getHistory(i).setXLabel("Frequency (w)");
		}
		this.transform = transform;
		this.nMeters = nMeters;
		current = new double[nMeters];	
	}

	public void updateValues() {
		current = transform.getData();
		for (int i=0;i<nMeters;i++) {
			currentValues[i] = current[i];
		}
	}

	public static double round(double x) {
			x=Math.round(x*1000);
			return x/1000;
		}

	public double[] getData() {
		return current;
	}

	public Dimension getDimension() {
		return Dimension.TIME;
	}

}

