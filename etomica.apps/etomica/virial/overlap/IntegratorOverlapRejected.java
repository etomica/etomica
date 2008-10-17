package etomica.virial.overlap;

import etomica.api.IRandom;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.util.Debug;
import etomica.virial.MeterVirialRejected;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlapRejected extends IntegratorManagerMC {

    public IntegratorOverlapRejected(IRandom random, IntegratorBox[] aIntegrators) {
        super(random);
        setNumSubSteps(1000);
        for (int i=0; i<aIntegrators.length; i++) {
            addIntegrator(aIntegrators[i]);
        }
        stepFreq = new double[nIntegrators];
        totNumSubSteps = new int[nIntegrators];
        setAdjustStepFreq(true);
        //there are no global moves
        setGlobalMoveInterval(Double.POSITIVE_INFINITY);

        // and hope nobody calls add/remove Integrators
    }
    
    /**
     * Sets the DataSource that retrieves data from both boxs and provides
     * information to the integrator about their progress.
     */
    public void setDSVO(DataSourceVirialOverlapRejected dataSource) {
        dsvo = dataSource;
        accumulators = dsvo.getAccumulators();
        if (integrators.length != accumulators.length) {
            throw new IllegalArgumentException("Must have the same number of integrators as accumulators\n");
        }
        for (int i=0; i<nIntegrators; i++) {
            if (doAdjustStepFreq) {
                stepFreq[i] = 1.0/nIntegrators;
            }
            totNumSubSteps[i] = 0;
        }
        minDiffLoc = -1;
    }
    
    public void reset() throws ConfigurationOverlapException {
        setDSVO(dsvo);
        super.reset();
    }

    /**
     * Sets the number of total number of integrator steps the sub-integrators
     * perform for every step of the overlap integrator.  Default value is 1000.
     */
    public void setNumSubSteps(int n) {
        numSubSteps = n;
    }
    
    /**
     * Sets whether to adjust the number of relative number of steps for each
     * sub-integrator.  Default is true.
     */
    public void setAdjustStepFreq(boolean b) {
        doAdjustStepFreq = b;
    }

    /**
     * Sets whether to adjust the number of relative number of steps for each
     * sub-integrator.  Default is true.
     */
    public boolean isAdjustStepFreq() {
        return doAdjustStepFreq;
    }

    // Override superclass so we can run the integrators for different lenghts
    // of time.  There are no global moves.
    public void doStepInternal() {
        for (int i=0; i<nIntegrators; i++) {
            int iSubSteps = numSubSteps/100 + (int)(numSubSteps*(1-0.01*nIntegrators) * stepFreq[i]); 
            for (int j=0; j<iSubSteps; j++) {
                integrators[i].doStep();
            }
            totNumSubSteps[i] += iSubSteps;
        }
        if (Debug.ON && Debug.DEBUG_NOW) {
            MeterVirialRejected[] meters = dsvo.getMeters();
            int newMinDiffLoc = dsvo.minDiffLocation();
            int nBennetPoints = meters[0].getNBennetPoints();
            if (nBennetPoints>1) System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+meters[0].getBennetBias(newMinDiffLoc)
                    +" ("+(((DataGroup)accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index).getValue(newMinDiffLoc+1)/
                            ((DataGroup)accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index).getValue(newMinDiffLoc+1))+")");
            System.out.print("Bennet bias ");
            for (int j=0; j<meters[0].getNBennetPoints(); j++) {
                System.out.print(meters[0].getBennetBias(j)+" ");
            }
            System.out.print("\n");
            for (int i=0; i<nIntegrators; i++) {
//                System.out.print("rel. uncertainty "+i+" ");
//                for (int j=0; j<accumulators[i].getNBennetPoints(); j++) {
//                    DataGroup data = (DataGroup)accumulators[i].getData(j);
//                    System.out.print(((DataDoubleArray)data.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]
//                                           /((DataDoubleArray)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]+" ");
//                }
                System.out.print("Bennet avg "+i+" ");
                for (int j=0; j<meters[i].getNBennetPoints(); j++) {
                    System.out.print(((DataGroup)accumulators[i].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index).getValue(j+1)+" ");
                }
                System.out.print("\n");
            }

        }
        if (doAdjustStepFreq) {
            adjustStepFreq();
        }
    }
    
    /**
     * Adjusts the frequency with which steps for each sub-integrator are
     * performed.  Integrators handling clusters with high errors are favored over
     * ones with low errors.
     */
    public void adjustStepFreq() {
        int newMinDiffLoc = dsvo.minDiffLocation();
        MeterVirialRejected[] meters = dsvo.getMeters();
        int nBennetPoints = meters[0].getNBennetPoints();
        if (newMinDiffLoc != minDiffLoc && nBennetPoints>1) {
            System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+meters[0].getBennetBias(newMinDiffLoc)
                    +" ("+(((DataGroup)accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index).getValue(newMinDiffLoc+1)/
                            ((DataGroup)accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index).getValue(newMinDiffLoc+1))+")");
        }
        minDiffLoc = newMinDiffLoc;

        double refError = ((DataGroup)accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index).getValue(minDiffLoc+1);
        double refErrorRatio = refError/Math.abs(((DataGroup)accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(minDiffLoc+1));
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println(0+" "+Math.abs(refError)+" "+Math.abs(refError/refErrorRatio));
        }
        if (Double.isNaN(refErrorRatio) || refErrorRatio > 1) {
            refErrorRatio = 1;
        }

        double targetError = ((DataGroup)accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index).getValue(minDiffLoc+1);
        double targetErrorRatio = targetError/Math.abs(((DataGroup)accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index).getValue(minDiffLoc+1));
        if (Double.isNaN(targetErrorRatio) || targetErrorRatio > 1) {
            targetErrorRatio = 1;
        }

        stepFreq[0] = 1.0 / (1 + targetErrorRatio/refErrorRatio * Math.sqrt(((double)totNumSubSteps[1]) / totNumSubSteps[0]));
        stepFreq[1] = 1.0 - stepFreq[0];

        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.print("freq ");
            for (int i=0; i<nIntegrators; i++) {
                System.out.print(stepFreq[i]+" ");
            }
            System.out.print("\n");
            System.out.print("steps ");
            for (int i=0; i<nIntegrators; i++) {
                System.out.print((double)totNumSubSteps[i]/(totNumSubSteps[0]+totNumSubSteps[1])+" ");
            }
            System.out.print("\n");
        }
    }
    
    public void setStepFreq0(double freq) {
        stepFreq[0] = freq;
        stepFreq[1] = 1-freq;
    }
    
    public double getStepFreq0() {
        return stepFreq[0];
    }
    
    public double getActualStepFreq0() {
        return totNumSubSteps[0]/(double)(totNumSubSteps[0]+totNumSubSteps[1]);
    }
    
    private static final long serialVersionUID = 1L;
    private final double[] stepFreq;
    private int numSubSteps;
    private int[] totNumSubSteps;
    private int minDiffLoc;
    private AccumulatorRatioAverage[] accumulators;
    private DataSourceVirialOverlapRejected dsvo;
    private boolean doAdjustStepFreq;
}
