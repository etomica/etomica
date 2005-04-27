package etomica.virial.overlap;

import etomica.Atom;
import etomica.Debug;
import etomica.Integrator;
import etomica.IntegratorIntervalEvent;
import etomica.PotentialMaster;
import etomica.data.AccumulatorAverage;
import etomica.utility.Arrays;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends Integrator {

    public IntegratorOverlap(PotentialMaster potentialMaster, Integrator[] aIntegrators, AccumulatorVirialOverlapSingleAverage[] virialAccumulators) {
        super(potentialMaster);
        numIntegrators = aIntegrators.length;
        integrators = aIntegrators;
        setNumSubSteps(1000);
        stepFreq = new double[numIntegrators];
        accumulators = virialAccumulators;
        errors = new double[numIntegrators];
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
            int iSubSteps = 10 + (int)(numSubSteps * stepFreq[i]); 
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
                    double[][] data = (double[][])accumulators[i].getTranslator().fromArray(accumulators[i].getData(j));
                    System.out.print(data[AccumulatorAverage.ERROR.index][1]/data[AccumulatorAverage.AVERAGE.index][1]+" ");
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
        if (newMinDiffLoc != minDiffLoc && nBennetPoints>1) System.out.println("target minDiffLoc = "+newMinDiffLoc+" refPref "+accumulators[0].getBennetBias(nBennetPoints-newMinDiffLoc-1));
        minDiffLoc = newMinDiffLoc;
        for (int i=0; i<numIntegrators; i++) {
            double[][] data;
            if (i==1) {
                data = (double[][])accumulators[i].getTranslator().fromArray(accumulators[i].getData(minDiffLoc));
            }
            else {
                data = (double[][])accumulators[i].getTranslator().fromArray(accumulators[i].getData(nBennetPoints-minDiffLoc-1));
            }
            errors[i] = data[AccumulatorAverage.ERROR.index][1] / data[AccumulatorAverage.AVERAGE.index][1];
//            System.out.println(i+" errors "+errors[i]);
//            System.out.print("Bennet "+i+" ");
//            for (int j=0; j<accumulators[i].getNBennetPoints(); j++) {
//                double[][] allData = (double[][])accumulators[i].getTranslator().fromArray(accumulators[i].getData(j));
//                System.out.print(allData[AccumulatorAverage.ERROR.index][1]/allData[AccumulatorAverage.AVERAGE.index][1]+" ");
//            }
//            System.out.print("\n");
            if (!Double.isNaN(errors[i]) && errors[i] > 0.0) {
                errors[i] *= errors[i];
            }
            else {
                errors[i] = 1.0;
            }
            stepFreq[i] = errors[i]*totNumSubSteps[i] / errors[0];
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
        System.out.print("freq ");
        for (int i=0; i<numIntegrators; i++) {
            stepFreq[i] /= freqTotal;
            System.out.print(stepFreq[i]+" ");
        }
        System.out.print("\n");
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

    public void reset() {
        for (int i=0; i<numIntegrators; i++) {
            integrators[i].reset();
        }
    }

    public Object makeAgent(Atom a) {
        return null;
    }

    private int numIntegrators;
    private Integrator[] integrators;
    private double[] stepFreq;
    private int numSubSteps;
    private int[] totNumSubSteps;
    private int minDiffLoc=-1;
    private AccumulatorVirialOverlapSingleAverage[] accumulators;
    private double[] errors;
    private IntegratorIntervalEvent[] intervalEvent;
    private DataSourceVirialOverlap dsvo;
    private boolean doAdjustStepFreq;
}
