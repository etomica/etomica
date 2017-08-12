/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.overlap;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorManagerMC;
import etomica.util.Debug;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends IntegratorManagerMC {

    public IntegratorOverlap(Integrator[] aIntegrators) {
        super(null);
        setNumSubSteps(1000);
        for (int i=0; i<aIntegrators.length; i++) {
            addIntegrator(aIntegrators[i]);
        }
        stepFreq = new double[nIntegrators];
        totNumSubSteps = new long[nIntegrators];
        setAdjustStepFreq(true);
        //there are no global moves
        setGlobalMoveInterval(Double.POSITIVE_INFINITY);

        // and hope nobody calls add/remove Integrators
    }
    
    /**
     * Sets the DataSource that retrieves data from both boxs and provides
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
    
    public void reset() {
        super.reset();
        if (dsvo != null) {
            setDSVO(dsvo);
        }
    }

    /**
     * Sets the number of total number of integrator steps the sub-integrators
     * perform for every step of the overlap integrator.  Default value is 1000.
     */
    public void setNumSubSteps(int n) {
        numSubSteps = n;
    }
    
    /**
     * Retruns the number of total number of integrator steps the sub-integrators
     * perform for every step of the overlap integrator.  Default value is 1000.
     */
    public long getNumSubSteps() {
        return numSubSteps;
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
    protected void doStepInternal() {
        for (int i=0; i<nIntegrators; i++) {
            long iSubSteps = (int)(numSubSteps*stepFreq[i]);
            if (doAdjustStepFreq) {
                // if we're internally adjusting the step fractions, require at least 1%
                iSubSteps = numSubSteps/100 + (int)(numSubSteps*(1-0.01*nIntegrators) * stepFreq[i]);
            }
            for (long j=0; j<iSubSteps; j++) {
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
//        int nBennetPoints = accumulators[0].getNBennetPoints();
//        if (newMinDiffLoc != minDiffLoc && nBennetPoints>1) {
//            System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc)
//                    +" ("+(accumulators[0].getBennetAverage(newMinDiffLoc)/accumulators[1].getBennetAverage(newMinDiffLoc))+")");
//        }
        minDiffLoc = newMinDiffLoc;

        DataGroup data = (DataGroup)accumulators[0].getData(minDiffLoc);
        double refError = ((DataDoubleArray)data.getData(accumulators[0].RATIO_ERROR.index)).getData()[1];
        double refErrorRatio = refError/Math.abs(((DataDoubleArray)data.getData(accumulators[0].RATIO.index)).getData()[1]);
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("0 "+Math.abs(refError)+" "+Math.abs(refError/refErrorRatio));
        }
        if (Double.isNaN(refErrorRatio) || refErrorRatio > 1) {
            // if we don't have enough data to calc the error, assume it's 100%
            // if apparent error is > 100%, cap it there (>100% just means our estimate for ratio is bad)
            refErrorRatio = 1.0;
        }

        data = (DataGroup)accumulators[1].getData(minDiffLoc);
        double targetError = ((DataDoubleArray)data.getData(accumulators[0].RATIO_ERROR.index)).getData()[1];
        double targetErrorRatio = targetError/Math.abs(((DataDoubleArray)data.getData(accumulators[0].RATIO.index)).getData()[1]);
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("1 "+Math.abs(targetError)+" "+Math.abs(targetError/targetErrorRatio));
        }
        if (Double.isNaN(targetErrorRatio) || targetErrorRatio > 1) {
            targetErrorRatio = 1.0;
        }
        stepFreq[0] = 1.0 / (1 + targetErrorRatio/refErrorRatio * Math.sqrt(((double)totNumSubSteps[1]) / totNumSubSteps[0]));
        stepFreq[1] = 1.0 - stepFreq[0];

        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("error ratios "+refErrorRatio+" "+targetErrorRatio);
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
    private long numSubSteps;
    private long[] totNumSubSteps;
    private int minDiffLoc;
    private AccumulatorVirialOverlapSingleAverage[] accumulators;
    private DataSourceVirialOverlap dsvo;
    private boolean doAdjustStepFreq;
}
