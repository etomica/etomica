package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataSourceScalar;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.units.Dimension;

/**
 * Measures ratio of two cluster integrals using overlap sampling.  The resulting ratio is 
 * formed from the ratio of a target and reference ratio from two different sub-simulations. 
 */
public class DataSourceVirialOverlap extends DataSourceScalar {

	private AccumulatorVirialOverlapSingleAverage refAccumulator, targetAccumulator;
	private final int nBennetPoints;
	
	public DataSourceVirialOverlap(AccumulatorVirialOverlapSingleAverage aRefAccumulator, 
			AccumulatorVirialOverlapSingleAverage aTargetAccumulator) {
        super("Virial Overlap Ratio",Dimension.FRACTION);
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
		nBennetPoints = aRefAccumulator.getNBennetPoints();
	}

    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the optimal value of the Bennet
     * parameter.
     */
	public double getDataAsScalar() {
		return getAverage(minDiffLocation());
	}
	
    public int getDataLength() {
        return 1;
    }
    
    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the given value of the Bennet
     * parameter.
     */
	public double getAverage(int iParam) {
        double targetAvg = ((DataDoubleArray)((DataGroup)targetAccumulator.getData(iParam)).getData(AccumulatorRatioAverage.RATIO.index)).getData()[1];
        double refAvg = ((DataDoubleArray)((DataGroup)refAccumulator.getData(nBennetPoints-iParam-1)).getData(AccumulatorRatioAverage.RATIO.index)).getData()[1];
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
                /targetAccumulator.getBennetAverage(0)-1);
		for (int i=1; i<nBennetPoints; i++) {
            // use average instead of sum to avoid race condition with overlap integrator
            double newDiff = Math.abs(refAccumulator.getBennetAverage(nBennetPoints-i-1)/
                             targetAccumulator.getBennetAverage(i)-1);
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
		double avg = getAverage(iParam);
        DataGroup dataGroup = (DataGroup)refAccumulator.getData(nBennetPoints-iParam-1);
		double refErr = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1];
        double refAvg = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.RATIO.index)).getData()[1];
        dataGroup = (DataGroup)targetAccumulator.getData(iParam);
		double targetErr = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1];
		double targetAvg = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.RATIO.index)).getData()[1];
		return Math.abs(avg)*Math.sqrt(refErr*refErr/(refAvg*refAvg)+targetErr*targetErr/(targetAvg*targetAvg));
	}

}