package etomica.virial.overlap;

import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.IntegratorPhase;
import etomica.potential.PotentialMaster;
import etomica.util.Debug;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends IntegratorManagerMC {

    public IntegratorOverlap(PotentialMaster potentialMaster, IntegratorPhase[] aIntegrators, AccumulatorVirialOverlapSingleAverage[] virialAccumulators) {
        super(potentialMaster);
        setNumSubSteps(1000);
        for (int i=0; i<aIntegrators.length; i++) {
            addIntegrator(aIntegrators[i]);
        }
        stepFreq = new double[nIntegrators];
        accumulators = virialAccumulators;
        totNumSubSteps = new int[nIntegrators];
        setAdjustStepFreq(true);
        //there are no global moves
        setGlobalMoveInterval(Double.POSITIVE_INFINITY);

        if (aIntegrators.length != accumulators.length) {
            throw new IllegalArgumentException("Must have the same number of integrators as accumulators\n");
        }
        // and hope nobody calls add/remove Integrators
    }
    
    /**
     * Sets the DataSource that retrieves data from both phases and provides
     * information to the integrator about their progress.
     */
    public void setDSVO(DataSourceVirialOverlap dataSource) {
        dsvo = dataSource;
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

    // Override superclass so we can run the integrators for different lenghts
    // of time.  There are no global moves.
    public void doStep() {
        for (int i=0; i<nIntegrators; i++) {
            int iSubSteps = 10 + (int)((numSubSteps-10*nIntegrators) * stepFreq[i]); 
            for (int j=0; j<iSubSteps; j++) {
                integrators[i].doStep();
                integrators[i].fireIntervalEvent(intervalEvents[i]);
            }
            totNumSubSteps[i] += iSubSteps;
        }
        if (Debug.ON && Debug.DEBUG_NOW) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            int nBennetPoints = accumulators[0].getNBennetPoints();
            if (nBennetPoints>1) System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc));
            for (int i=0; i<nIntegrators; i++) {
                System.out.print("Bennet "+i+" ");
                for (int j=0; j<accumulators[i].getNBennetPoints(); j++) {
                    DataGroup data = (DataGroup)accumulators[i].getData(j);
                    System.out.print(((DataDoubleArray)data.getData(AccumulatorAverage.ERROR.index)).getData()[1]
                                           /((DataDoubleArray)data.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]+" ");
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
        double freqTotal = 0;
        int newMinDiffLoc = dsvo.minDiffLocation();
        int nBennetPoints = accumulators[0].getNBennetPoints();
        if (newMinDiffLoc != minDiffLoc && nBennetPoints>1) {
            System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc));
        }
        minDiffLoc = newMinDiffLoc;
        for (int i=0; i<nIntegrators; i++) {
            DataGroup data;
            data = (DataGroup)accumulators[i].getData(minDiffLoc);
            double error = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1];
            data = (DataGroup)accumulators[1-i].getData(minDiffLoc);
            double otherRatio = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.RATIO.index)).getData()[1];
//            System.out.println(i+" errors "+errors[i]);
//            System.out.print("Bennet "+i+" ");
//            for (int j=0; j<accumulators[i].getNBennetPoints(); j++) {
//                double[][] allData = (double[][])accumulators[i].getTranslator().fromArray(accumulators[i].getData(j));
//                System.out.print(allData[AccumulatorAverage.ERROR.index][1]/allData[AccumulatorAverage.AVERAGE.index][1]+" ");
//            }
//            System.out.print("\n");
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println(i+" "+Math.abs(error)+" "+Math.abs(otherRatio));
            }
            if (!Double.isNaN(error) && !Double.isNaN(otherRatio)) {
                error *= error;
                otherRatio *= otherRatio;
            }
            else {
                error = 1.0;
                otherRatio = 1.0;
            }
            stepFreq[i] = otherRatio*error*totNumSubSteps[i];
            freqTotal += stepFreq[i];
        }
//        System.out.print("ratio error ");
//        for (int j=0; j<accumulators[0].getNBennetPoints(); j++) {
//            System.out.print(Math.sqrt(allErrors[nBennetPoints-j-1][0]*allErrors[nBennetPoints-j-1][0]+allErrors[j][1]*allErrors[j][1])+" ");
//        }
//        System.out.print("\n");
//        System.out.println("avgs "+accumulators[0].getData(AccumulatorVirialAverage.AVERAGE_RATIO)[1]+" "+accumulators[1].getData(AccumulatorVirialAverage.AVERAGE_RATIO)[1]);
//        System.out.println("ratio "+accumulators[0].getData(AccumulatorVirialAverage.AVERAGE_RATIO)[1]/accumulators[1].getData(AccumulatorVirialAverage.AVERAGE_RATIO)[1]);
//        System.out.println("errors "+errors[0]+" "+errors[1]);
        // normalize
        for (int i=0; i<nIntegrators; i++) {
            stepFreq[i] /= freqTotal;
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
    
    private final double[] stepFreq;
    private int numSubSteps;
    private int[] totNumSubSteps;
    private int minDiffLoc=-1;
    private AccumulatorVirialOverlapSingleAverage[] accumulators;
    private DataSourceVirialOverlap dsvo;
    private boolean doAdjustStepFreq;
}
