/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterPressureHardTensor;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.meter.MeterStructureFactor;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.util.ParseArgs;

public class GlassClustering {
    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 12;
            params.nB = 12;
            params.density = 1.6;
            params.D = 3;
            params.temperature = 1.0;
            params.numSteps = 1000000;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);
        System.out.println(params.D + "D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        System.out.println("T = " + params.temperature);
        System.out.println(params.numSteps + " MD steps after " + params.numSteps / 10 + " equilibaration steps");

        //Equilibration
        double temperature0 = Math.max(params.temperatureMelt, params.temperature);
        if (temperature0 > params.temperature) System.out.println("Equilibrating at T=" + temperature0);
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(temperature0);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps / 10);


        if (temperature0 > params.temperature) {
            System.out.println("Equilibrating at T=" + params.temperature);
            sim.integrator.setTemperature(params.temperature);
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps / 10);

        }

        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);
IDataSource pTensorMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pTensorMeter = new MeterPressureTensorFromIntegrator(sim.getSpace());
            ((MeterPressureTensorFromIntegrator) pTensorMeter).setIntegrator((IntegratorVelocityVerlet) sim.integrator);
        } else {
            pTensorMeter = new MeterPressureHardTensor(sim.getSpace());
            ((MeterPressureHardTensor) pTensorMeter).setIntegrator((IntegratorHard) sim.integrator);
            new MeterPressureHard((IntegratorHard) sim.integrator);
        }

        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        DataPumpListener pTensorPump = new DataPumpListener(pTensorMeter, tracer, params.interval);

        sim.integrator.getEventManager().addListener(pTensorPump);

        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);

        double cut = params.cut;
        System.out.println("using Sfac cutoff " + cut);

        MeterStructureFactor meterSFacCluster = new MeterStructureFactor(sim.box, cut);
        DataClusterWriter clusterWriter = new DataClusterWriter(sim.box, "sfac.bin");
        DataPumpListener pumpSFacCluster = new DataPumpListener(meterSFacCluster, clusterWriter, params.interval);
        sim.integrator.getEventManager().addListener(pumpSFacCluster);
        pFork.addDataSink(clusterWriter.makePressureSink());
        double vB = sim.getSpace().powerD(sim.sigmaB);
        ((MeterStructureFactor.AtomSignalSourceByType) meterSFacCluster.getSignalSource()).setAtomTypeFactor(sim.speciesB.getAtomType(0), vB);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(params.numSteps / 100 / params.interval);
        pFork.addDataSink(accP);

        //Run
        double time0 = System.nanoTime();
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps);

        //P
        double time1 = System.nanoTime();

        clusterWriter.closeFile();

        IData pData = accP.getData();
        double pAvg = pData.getValue(accP.AVERAGE.index);
        double pErr = pData.getValue(accP.ERROR.index);
        double pCor = pData.getValue(accP.BLOCK_CORRELATION.index);

        System.out.println("P: " + pAvg + " err: " + pErr + " cor: " + pCor);

        System.out.println("\ntime: " + (time1 - time0) / 1e9);
    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numSteps = 1000000;
        public double temperatureMelt = 0;
        public int interval = 10;
        public double cut = 8;
    }
}