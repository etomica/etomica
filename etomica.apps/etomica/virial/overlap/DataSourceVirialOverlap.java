package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataSourceAdapter;
import etomica.units.Dimension;

/**
 * Measures ratio of two cluster integrals using overlap sampling.  The resulting ratio is 
 * formed from the ratio of a target and reference ratio from two different sub-simulations. 
 */
public class DataSourceVirialOverlap extends DataSourceAdapter {

	private AccumulatorVirialOverlapSingleAverage refAccumulator, targetAccumulator;
	private final int nBennetPoints;
	
	public DataSourceVirialOverlap(AccumulatorVirialOverlapSingleAverage aRefAccumulator, 
			AccumulatorVirialOverlapSingleAverage aTargetAccumulator) {
        super(Dimension.FRACTION);
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
		nBennetPoints = aRefAccumulator.getNBennetPoints();
	}

    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the optimal value of the Bennet
     * parameter.
     */
	public double[] getData() {
		return new double[] {getData(minDiffLocation())};
	}
	
    public int getDataLength() {
        return 1;
    }
    
    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the given value of the Bennet
     * parameter.
     */
	public double getData(int iParam) {
        double targetAvg = ((double[][])targetAccumulator.getTranslator().fromArray(targetAccumulator.getData(iParam)))[AccumulatorRatioAverage.RATIO.index][1];
        double refAvg = ((double[][])refAccumulator.getTranslator().fromArray(refAccumulator.getData(nBennetPoints-iParam-1)))[AccumulatorRatioAverage.RATIO.index][1];
        return refAvg/targetAvg;
	}
	
	/**
     * Returns the index of the Bennet parameter which minimizes the differences 
     * between the Bennet "sums" from the target and references accumulators.  This
     * parameter should be optimal for overlap sampling.
	 */
    public int minDiffLocation() {
		int minDiffLoc = 0;
		double minDiff = Math.abs(refAccumulator.getBennetAverage(nBennetPoints-1)
                -targetAccumulator.getBennetAverage(0));
		for (int i=1; i<nBennetPoints; i++) {
            // use average instead of sum to avoid race condition with overlap integrator
            double newDiff = Math.abs(refAccumulator.getBennetAverage(nBennetPoints-i-1)
                             -targetAccumulator.getBennetAverage(i));
			if (newDiff < minDiff) {
				minDiffLoc = i;
				minDiff = newDiff;
			}
		}
		return minDiffLoc;
	}
	
    /**
     * Returns the error in the ratio of the reference to target 
     * overlap-to-virial ratios (which reduces to target/reference) 
     * for the optimal value of the Bennet parameter.
     */
	public double getError() {
		return getError(minDiffLocation());
	}

    /**
     * Returns the error in the ratio of the reference to target 
     * overlap-to-virial ratios (which reduces to target/reference) 
     * for the given value of the Bennet parameter.
     */
	public double getError(int iParam) {
		double avg = getData(iParam);
        double[][] data = (double[][])refAccumulator.getTranslator().fromArray(refAccumulator.getData(nBennetPoints-iParam-1));
		double refErr = data[AccumulatorRatioAverage.RATIO_ERROR.index][1];
        double refAvg = data[AccumulatorRatioAverage.RATIO.index][1];
        data = (double[][])targetAccumulator.getTranslator().fromArray(targetAccumulator.getData(iParam));
		double targetErr = data[AccumulatorRatioAverage.RATIO_ERROR.index][1];
		double targetAvg = data[AccumulatorRatioAverage.RATIO.index][1];
		return Math.abs(avg)*Math.sqrt(refErr*refErr/(refAvg*refAvg)+targetErr*targetErr/(targetAvg*targetAvg));
	}

}