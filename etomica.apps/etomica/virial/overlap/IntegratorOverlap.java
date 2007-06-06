package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.IntegratorPhase;
import etomica.util.Debug;
import etomica.util.IRandom;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends IntegratorManagerMC {

    public IntegratorOverlap(IRandom random, IntegratorPhase[] aIntegrators) {
        super(null, random);
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
     * Sets the DataSource that retrieves data from both phases and provides
     * information to the integrator about their progress.
     */
    public void setDSVO(DataSourceVirialOverlap dataSource) {
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
            int newMinDiffLoc = dsvo.minDiffLocation();
            int nBennetPoints = accumulators[0].getNBennetPoints();
            if (nBennetPoints>1) System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc)
                    +" ("+(accumulators[0].getBennetAverage(newMinDiffLoc)/accumulators[1].getBennetAverage(newMinDiffLoc))+")");
            System.out.print("Bennet bias ");
            for (int j=0; j<accumulators[0].getNBennetPoints(); j++) {
                System.out.print(accumulators[0].getBennetBias(j)+" ");
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
                for (int j=0; j<accumulators[i].getNBennetPoints(); j++) {
                    System.out.print(accumulators[i].getBennetAverage(j)+" ");
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
        int nBennetPoints = accumulators[0].getNBennetPoints();
        if (newMinDiffLoc != minDiffLoc && nBennetPoints>1) {
            System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc)
                    +" ("+(accumulators[0].getBennetAverage(newMinDiffLoc)/accumulators[1].getBennetAverage(newMinDiffLoc))+")");
        }
        minDiffLoc = newMinDiffLoc;
        for (int i=0; i<nIntegrators; i++) {
            DataGroup data;
            data = (DataGroup)accumulators[i].getData(minDiffLoc);
            double error = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1];
            double ratio = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println(i+" "+Math.abs(error)+" "+Math.abs(ratio));
            }
            if (Double.isNaN(error) || Double.isNaN(ratio)) {
                error = 1.0;
                ratio = 1.0;
            }
            double sum = 1.0;
            for (int j=0; j<nIntegrators; j++) {
                if (j==i) continue;
                data = (DataGroup)accumulators[j].getData(minDiffLoc);
                double jRatio = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
                double jError = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1];
                if (Double.isNaN(jError) || Double.isNaN(jRatio)) {
                    jError = 1.0;
                    jRatio = 1.0;
                }
                double r = (ratio/error)*(jError/jRatio);
                sum += r * r * totNumSubSteps[j] / totNumSubSteps[i];
            }
            stepFreq[i] = 1.0 / sum;
        }
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
    private AccumulatorVirialOverlapSingleAverage[] accumulators;
    private DataSourceVirialOverlap dsvo;
    private boolean doAdjustStepFreq;
}
