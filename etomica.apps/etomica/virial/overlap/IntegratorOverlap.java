package etomica.virial.overlap;

import etomica.atom.Atom;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.potential.PotentialMaster;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends Integrator {

    public IntegratorOverlap(PotentialMaster potentialMaster, Integrator[] aIntegrators, AccumulatorVirialOverlapSingleAverage[] virialAccumulators) {
        super(potentialMaster,0);
        numIntegrators = aIntegrators.length;
        integrators = aIntegrators;
        setNumSubSteps(1000);
        stepFreq = new double[numIntegrators];
        accumulators = virialAccumulators;
        intervalEvent = new IntegratorIntervalEvent[numIntegrators];
        totNumSubSteps = new int[numIntegrators];
        for (int i=0; i<numIntegrators; i++) {
            intervalEvent[i] = new IntegratorIntervalEvent(integrators[i],1);
        }
        setAdjustStepFreq(true);
        
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
    
    public void doStep() {
        for (int i=0; i<numIntegrators; i++) {
            int iSubSteps = 10 + (int)((numSubSteps-10*numIntegrators) * stepFreq[i]); 
            for (int j=0; j<iSubSteps; j++) {
                integrators[i].doStep();
                integrators[i].fireIntervalEvent(intervalEvent[i]);
            }
            totNumSubSteps[i] += iSubSteps;
        }
        if (Debug.ON && Debug.DEBUG_NOW) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            int nBennetPoints = accumulators[0].getNBennetPoints();
            if (nBennetPoints>1) System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(newMinDiffLoc));
            for (int i=0; i<numIntegrators; i++) {
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
            System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(nBennetPoints-newMinDiffLoc-1));
        }
        minDiffLoc = newMinDiffLoc;
        for (int i=0; i<numIntegrators; i++) {
            DataGroup data;
            if (i==1) {
                data = (DataGroup)accumulators[i].getData(minDiffLoc);
            }
            else {
                data = (DataGroup)accumulators[i].getData(nBennetPoints-minDiffLoc-1);
            }
            double error = ((DataDoubleArray)data.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1];
            if (i==1) {
                data = (DataGroup)accumulators[1-i].getData(nBennetPoints-minDiffLoc-1);
            }
            else {
                data = (DataGroup)accumulators[1-i].getData(minDiffLoc);
            }
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
        for (int i=0; i<numIntegrators; i++) {
            stepFreq[i] /= freqTotal;
        }
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.print("freq ");
            for (int i=0; i<numIntegrators; i++) {
                System.out.print(stepFreq[i]+" ");
            }
            System.out.print("\n");
            System.out.print("steps ");
            for (int i=0; i<numIntegrators; i++) {
                System.out.print((double)totNumSubSteps[i]/(totNumSubSteps[0]+totNumSubSteps[1])+" ");
            }
            System.out.print("\n");
        }
    }
    
    /**
     * Adds a sub-integrator to the overlap integrator.  Steps for each
     * integrator are performed every step of the overlap integrator.
     */
    public void addIntegrator(Integrator integrator) {
        integrators = (Integrator[])Arrays.addObject(integrators,integrator);
        numIntegrators++;
    }
    
    /**
     * Removes the given integrator from the list of integrators.  Returns
     * false if the given integrattor was not handled by this overlap integrator.
     */
    public boolean removeIntegrator(Integrator integrator) {
        integrators = (Integrator[])Arrays.removeObject(integrators,integrator);
        if (numIntegrators == integrators.length) {
            return false;
        }
        numIntegrators--;
        return true;
    }

    public void reset() throws ConfigurationOverlapException {
        ConfigurationOverlapException overlapException = null;
	    for(int i=0; i<numIntegrators; i++) {
            try {
                integrators[i].reset();
            }
            catch (ConfigurationOverlapException e) {
                if (overlapException == null) {
                    overlapException = e;
                }
            }
        }
        if (overlapException != null) {
            throw overlapException;
        }
        // don't call super.reset(), it tries to calculate the potential energy
    }

    public Object makeAgent(Atom a) {
        return null;
    }
    
    public void setEquilibrating(boolean flag) {
        for (int i=0; i<numIntegrators; i++) {
            integrators[i].setEquilibrating(flag);
        }
    }        

    public void setStepFreq0(double freq) {
        stepFreq[0] = freq;
        stepFreq[1] = 1-freq;
    }
    
    public double getStepFreq0() {
        return stepFreq[0];
    }
    
    private int numIntegrators;
    private Integrator[] integrators;
    private final double[] stepFreq;
    private int numSubSteps;
    private int[] totNumSubSteps;
    private int minDiffLoc=-1;
    private AccumulatorVirialOverlapSingleAverage[] accumulators;
    private IntegratorIntervalEvent[] intervalEvent;
    private DataSourceVirialOverlap dsvo;
    private boolean doAdjustStepFreq;
}
