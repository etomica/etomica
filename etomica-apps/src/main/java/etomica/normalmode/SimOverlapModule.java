/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.ActivityIntegrate2;
import etomica.action.activity.Controller;
import etomica.action.controller.Controller2;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.data.AccumulatorAverageCollapsingLog;
import etomica.data.DataPumpListener;
import etomica.data.DataSplitter;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.IDataSink;
import etomica.integrator.IntegratorMC;
import etomica.overlap.DataOverlap;
import etomica.overlap.DataOverlap.DataSourceOverlapLogAvg;
import etomica.overlap.IntegratorOverlap;
import etomica.overlap.MeterOverlap;

/**
 * Overlap sampling "module" to calculate free energy difference.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapModule {

    public SimOverlapModule(Box boxReference, Box boxTarget,
                            Integrator integratorReference, Integrator integratorTarget,
                            IAPIPotential potentialReference, IAPIPotential potentialTarget,
                            double temperature) {
        integrators = new Integrator[2];
        integrators[0] = integratorReference;
        integrators[1] = integratorTarget;
        accumulatorPumps = new DataPumpListener[2];
        meters = new MeterOverlap[2];

        MeterAPIPotentialEnergy meterTarget = new MeterAPIPotentialEnergy(potentialTarget);
        meterTarget.setBox(boxTarget);
        MeterAPIPotentialEnergy meterTargetInReference = new MeterAPIPotentialEnergy(potentialTarget);
        meterTargetInReference.setBox(boxReference);
        MeterAPIPotentialEnergy meterReference = new MeterAPIPotentialEnergy(potentialReference);
        meterReference.setBox(boxReference);
        MeterAPIPotentialEnergy meterReferenceInTarget = new MeterAPIPotentialEnergy(potentialReference);
        meterReferenceInTarget.setBox(boxTarget);
        integratorOverlap = new IntegratorOverlap(integrators);
        meters[0] = new MeterOverlap(meterReference, meterTargetInReference, temperature, true);
        meters[1] = new MeterOverlap(meterTarget, meterReferenceInTarget, temperature, false);

        final IDataSinkFactory dataSinkFactory = new IDataSinkFactory() {
            public IDataSink makeDataSink(int i) {
                return new AccumulatorAverageCollapsingLog();
            }
        };

        DataSplitter splitterReference = new DataSplitter();
        splitterReference.setDataSinkFactory(dataSinkFactory);
        accumulatorPumps[0] = new DataPumpListener(meters[0], splitterReference);
        integratorReference.getEventManager().addListener(accumulatorPumps[0]);
        splitterReference.setDataSinkFactory(dataSinkFactory);
        DataSourceOverlapLogAvg overlapAvgA = new DataOverlap.DataSourceOverlapAvgCollapsingSplit(splitterReference);

        DataSplitter splitterTarget = new DataSplitter();
        splitterTarget.setDataSinkFactory(dataSinkFactory);
        accumulatorPumps[1] = new DataPumpListener(meters[1], splitterTarget);
        integratorTarget.getEventManager().addListener(accumulatorPumps[1]);
        splitterTarget.setDataSinkFactory(dataSinkFactory);
        DataSourceOverlapLogAvg overlapAvgB = new DataOverlap.DataSourceOverlapAvgCollapsingSplit(splitterTarget);

        dataOverlap = new DataOverlap(overlapAvgA, overlapAvgB, meters[0]);
        
        integratorOverlap.setReferenceFracSource(dataOverlap);
        
        setAlpha(1.0, 30);
        
        controller = new Controller2();
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

    public DataOverlap getDataOverlap() {
        return dataOverlap;
    }

    public double getAlphaCenter() {
        return meters[0].getAlphaCenter();
    }

    public void setAlpha(double refPrefCenter, double span) {
        meters[0].setAlphaRange(refPrefCenter, span);
        meters[1].setAlphaRange(refPrefCenter, span);
    }

    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        double refPref = -1.0;
        if (fileName != null) {
            try { 
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+refPref);
                double span = meters[0].getAlphaSpan();
                setAlpha(refPref, span);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            controller.runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps/2);
            System.out.println("target equilibration finished");

            setAlpha(1, 60);
            controller.runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps);
            // in theory we could care about equilibrating the reference, so we'll do that too
            // but, in practice, the reference needs no equilibration

            refPref = dataOverlap.getOverlapAverageAndError()[0];
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            
            setAlpha(refPref,5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
        }

    }
    
    public void equilibrate(long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        controller.runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps);
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        // this is just to cause the averages to reset
        setAlpha(meters[0].getAlphaCenter(), meters[0].getAlphaSpan());
    }

    private static final long serialVersionUID = 1L;
    protected IntegratorOverlap integratorOverlap;
    protected DataOverlap dataOverlap;
    protected Integrator[] integrators;
    protected DataPumpListener[] accumulatorPumps;
    protected MeterOverlap[] meters;
    protected Controller2 controller;
    protected int targetDataInterval, referenceDataInterval;
}
