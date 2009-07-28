package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.units.Fraction;
import etomica.virial.MeterVirialRejected;

/**
 * Measures ratio of two cluster integrals using overlap sampling.  The resulting ratio is 
 * formed from the ratio of a target and reference ratio from two different sub-simulations. 
 */
public class DataSourceVirialOverlapRejected extends DataSourceScalar {

    private static final long serialVersionUID = 1L;
    protected AccumulatorRatioAverage refAccumulator, targetAccumulator;
    protected MeterVirialRejected refMeter, targetMeter;
	private final int nBennetPoints;
	
	public DataSourceVirialOverlapRejected(AccumulatorRatioAverage aRefAccumulator, 
			AccumulatorRatioAverage aTargetAccumulator, MeterVirialRejected aRefMeter,
			MeterVirialRejected aTargetMeter) {
        super("Virial Overlap Ratio",Fraction.DIMENSION);
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
		refMeter = aRefMeter;
		targetMeter = aTargetMeter;
        nBennetPoints = aRefMeter.getNBennetPoints();
	}
    
    public AccumulatorRatioAverage[] getAccumulators() {
        return new AccumulatorRatioAverage[]{refAccumulator, targetAccumulator};
    }

    public MeterVirialRejected[] getMeters() {
        return new MeterVirialRejected[]{refMeter, targetMeter};
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
        double targetAvg = ((DataGroup)targetAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(iParam+1);
        double refAvg = ((DataGroup)refAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(iParam+1);
        return refAvg/targetAvg;
	}
	
	/**
     * Returns the index of the Bennet parameter which minimizes the differences 
     * between the Bennet "sums" from the target and references accumulators.  This
     * parameter should be optimal for overlap sampling.
	 */
    public int minDiffLocation() {
		int minDiffLoc = 0;
		IData avgRefData = ((DataGroup)refAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        IData avgTargetData = ((DataGroup)targetAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        double ratio = avgRefData.getValue(0+1)/avgTargetData.getValue(0+1);
        double bias = refMeter.getBennetBias(0);
		double minDiff = ratio/bias + bias/ratio - 2;
		for (int i=1; i<nBennetPoints; i++) {
            ratio = avgRefData.getValue(i+1)/avgTargetData.getValue(i+1);
            bias = refMeter.getBennetBias(i);
            double newDiff = ratio/bias + bias/ratio - 2;
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
		double refErr = ((DataGroup)refAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index).getValue(iParam+1);
        double refAvg = ((DataGroup)refAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(iParam+1);
        double refRelErr = refErr/refAvg;
		double targetErr = ((DataGroup)targetAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index).getValue(iParam+1);
		double targetAvg = ((DataGroup)targetAccumulator.getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(iParam+1);
        double targetRelErr = targetErr/targetAvg;
		return Math.abs(avg)*Math.sqrt(refRelErr*refRelErr+targetRelErr*targetRelErr);
	}

    /**
     * Convenience method that resets reference and target accumulators
     */
    public void reset() {
        targetAccumulator.reset();
        refAccumulator.reset();
    }
}