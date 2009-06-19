package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.data.DataPumpListener;
import etomica.data.IEtomicaDataSource;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * Overlap sampling "module" to calculate free energy difference.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapModule {

    public SimOverlapModule(IntegratorBox[] integrators, IEtomicaDataSource meterTarget,
            IEtomicaDataSource meterReference, IEtomicaDataSource meterTargetInReference, IEtomicaDataSource meterReferenceInTarget,
            double temperature) {

        this.integrators = integrators;
        accumulatorPumps = new DataPumpListener[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(integrators);
        meters[1] = new MeterOverlap(meterTarget, meterReferenceInTarget, temperature);
        meters[0] = new MeterOverlap(meterReference, meterTargetInReference, temperature);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        
        setRefPref(1.0, 30);
        
        activityIntegrate = new ActivityIntegrate(integratorOverlap);
        controller = new Controller();
        controller.addAction(activityIntegrate);
    }
    
    public void setTargetDataInterval(int newTargetDataInterval) {
        targetDataInterval = newTargetDataInterval;
        accumulatorPumps[1].setInterval(targetDataInterval);
    }
    
    public int getTargetDataInterval() {
        return targetDataInterval;
    }

    public void setReferenceDataInterval(int newReferenceDataInterval) {
        referenceDataInterval = newReferenceDataInterval;
        accumulatorPumps[0].setInterval(referenceDataInterval);
    }

    public int getReferenceDataInterval() {
        return referenceDataInterval;
    }

    public IntegratorOverlap getIntegratorOverlap() {
        return integratorOverlap;
    }

    public void setIntegratorOverlap(IntegratorOverlap integratorOverlap) {
        this.integratorOverlap = integratorOverlap;
    }

    public DataSourceVirialOverlap getDsvo() {
        return dsvo;
    }

    public ActivityIntegrate getActivityIntegrate() {
        return activityIntegrate;
    }

    public double getRefPref() {
        return refPref;
    }

    public AccumulatorVirialOverlapSingleAverage[] getAccumulators() {
        return accumulators;
    }

    public Controller getController() {
        return controller;
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (iBox == 0) {
            newAccumulator.setBlockSize(100);
        }
        else {
            newAccumulator.setBlockSize(100);
        }
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPumpListener(meters[iBox],newAccumulator);
            integrators[iBox].getEventManager().addListener(accumulatorPumps[iBox]);
            if (iBox == 1) {
                accumulatorPumps[iBox].setInterval(targetDataInterval);
            }
            else {
                accumulatorPumps[iBox].setInterval(referenceDataInterval);
            }
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setDSVO(dsvo);
        }
    }
    
    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setRefPref(newRefPref,1);
    }
    
    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try { 
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps/2);
            controller.actionPerformed();
            controller.reset();
            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1,60);
            activityIntegrate.setMaxSteps(initSteps);
            controller.actionPerformed();
            controller.reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setRefPref(refPref,5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
            controller.reset();
        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        controller.actionPerformed();
        controller.reset();
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
            setRefPref(refPref,1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dsvo.reset();
        }
    }

    private static final long serialVersionUID = 1L;
    protected IntegratorOverlap integratorOverlap;
    protected DataSourceVirialOverlap dsvo;
    protected IntegratorBox[] integrators;
    protected ActivityIntegrate activityIntegrate;
    protected double refPref;
    protected AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected DataPumpListener[] accumulatorPumps;
    protected IEtomicaDataSource[] meters;
    protected Controller controller;
    protected int targetDataInterval, referenceDataInterval;
}
