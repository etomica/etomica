/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterPressureHardTensor;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;

public class GlassProdF {
    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 125;
            params.nB = 125;
            params.density = 1.4;
            params.D = 3;
            params.temperature = 0.1;
            params.numStepsEq = 100000;
            params.numSteps =   1000000;
            params.minDrFilter = 0.4;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential);
        System.out.println(params.D +"D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        double volume = sim.box.getBoundary().volume();
        int numAtoms = params.nA + params.nB;
        double rho= numAtoms/volume;
        System.out.println( params.numSteps + " MD steps after " + params.numStepsEq + " equilibaration steps");

        //Equilibration
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(params.temperature);
        sim.activityIntegrate.setMaxSteps(params.numStepsEq);
        sim.getController().actionPerformed();
        sim.getController().reset();

        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);
        sim.activityIntegrate.setMaxSteps(params.numSteps);

        //P
        IDataSource pTensorMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pTensorMeter = new MeterPressureTensorFromIntegrator(sim.getSpace());
            ((MeterPressureTensorFromIntegrator) pTensorMeter).setIntegrator((IntegratorVelocityVerlet) sim.integrator);
        } else {
            pTensorMeter = new MeterPressureHardTensor(sim.getSpace());
            ((MeterPressureHardTensor) pTensorMeter).setIntegrator((IntegratorHard) sim.integrator);
            new MeterPressureHard((IntegratorHard) sim.integrator);
        }


        DataFork pTensorFork = new DataFork();
        DataPumpListener pTensorPump = new DataPumpListener(pTensorMeter, pTensorFork);
        AccumulatorAverageFixed pTensorAccumulator = new AccumulatorAverageFixed(1);
        pTensorFork.addDataSink(pTensorAccumulator);

        sim.integrator.getEventManager().addListener(pTensorPump);



        long blocksize = params.numSteps/100;
        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        pFork.addDataSink(pAccumulator);



        //MSD
        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSD.setEnabled(true);

        //Fs
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        configStorageMSD.addListener(meterFs);

        //F
        DataSourceFOld meterF = new DataSourceFOld(configStorageMSD);
        configStorageMSD.addListener(meterF);

        sim.integrator.getEventManager().addListener(configStorageMSD);

        //Run
        double time0 = System.currentTimeMillis();
        sim.getController().actionPerformed();

        int[] pIndex;
        if(params.D==2){
            pIndex = new int[] {1};
        }else{
            pIndex = new int[] {1,2,5};
        }

        //Pressure
        DataGroup dataP = (DataGroup)pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg  = dataPAvg.getValue(0);
        double pErr  = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);

        String filenameF, filenameFs;

        if(sim.potentialChoice == SimGlass.PotentialChoice.HS){
            double phi;
            if(params.D == 2){
                phi = Math.PI/4*(params.nA+params.nB/(1.4*1.4))/volume;
            }else{
                phi = Math.PI/6*(params.nA+params.nB/(1.4*1.4*1.4))/volume;
            }
            System.out.println("rho: " + params.density + "  phi: " + phi+"\n");
            System.out.println("Z: " + pAvg/params.density/params.temperature +"  "+ pErr/params.density/params.temperature  +"  cor: "+pCorr);
            filenameFs = String.format("fsRho%1.3f.out", rho);
            filenameF = String.format("fRho%1.3f.out", rho);
        }else{
            System.out.println("rho: " + params.density+"\n");
            filenameFs = String.format("fsRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameF = String.format("fRho%1.3fT%1.3f.out",  rho, params.temperature);
        }
        System.out.println("P: " + pAvg +"  "+ pErr +"  cor: "+pCorr);



        //MSD
        FileWriter fileWriterFs, fileWriterF;
        try {
            fileWriterFs  = new FileWriter(filenameFs, false);
            fileWriterF   = new FileWriter(filenameF,  false);
            DataDoubleArray x = meterF.getIndependentData(0);
            for (int i=0; i<meterF.getData().getLength(); i++){
                double xi = x.getValue(i);
                double yiF  = meterF.getData().getValue(i);
                double yiFs = meterFs.getData().getValue(i);
                if( !Double.isNaN(yiF) ){
                    fileWriterF.write(xi + " " + yiF + "\n");
                    fileWriterFs.write(xi + " " + yiFs + "\n");
                }
            }
            fileWriterF.close();
            fileWriterFs.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        double time1 = System.currentTimeMillis();
        System.out.println("\ntime: " + (time1-time0)/1000/3600 + " hrs");
    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;
    }
}