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
public class GlassProd {
    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 100;
            params.nB = 100;
            params.density = 1.;
            params.D = 3;
            params.temperature = 1.0;
            params.numStepsEq = 100000;
            params.numSteps =   1000000;
            params.log2StepS = 5;
            params.log2StepE = 20;
            params.minDrFilter = 0.5;
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

        //G
        IDataSource pTensorMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pTensorMeter = new MeterPressureTensorFromIntegrator(sim.getSpace());
            ((MeterPressureTensorFromIntegrator) pTensorMeter).setIntegrator((IntegratorVelocityVerlet) sim.integrator);
        } else {
            pTensorMeter = new MeterPressureHardTensor(sim.getSpace());
            ((MeterPressureHardTensor) pTensorMeter).setIntegrator((IntegratorHard) sim.integrator);
            new MeterPressureHard((IntegratorHard) sim.integrator);
        }

        //Viscosity
        DataFork pTensorFork = new DataFork();
        AccumulatorPTensor pTensorAccumVisc = new AccumulatorPTensor(sim.integrator, sim.integrator.getTimeStep());
        pTensorAccumVisc.setEnabled(true);
        pTensorFork.addDataSink(pTensorAccumVisc);
        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork);


        AccumulatorAverageFixed pTensorAccumulator = new AccumulatorAverageFixed(1);
        pTensorFork.addDataSink(pTensorAccumulator);

        sim.integrator.getEventManager().addListener(pTensorAccumViscPump);

        long blocksize = params.numSteps/100;
        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        pFork.addDataSink(pAccumulator);

        AccumulatorAverageFixed tAccumulator = null;
        if(sim.potentialChoice != SimGlass.PotentialChoice.HS){
            MeterTemperature tMeter = new MeterTemperature(sim, sim.box, params.D);
            tAccumulator = new AccumulatorAverageFixed(blocksize);
            DataPumpListener tPump = new DataPumpListener(tMeter, tAccumulator, 10);
            tAccumulator.addDataSink(pTensorAccumVisc.makeTemperatureSink(),  new AccumulatorAverage.StatType[]{tAccumulator.AVERAGE});
            sim.integrator.getEventManager().addListener(tPump);
        }


        //MSD
        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSD.setEnabled(true);
        DataSourceMSD meterMSD = new DataSourceMSD(configStorageMSD);
        configStorageMSD.addListener(meterMSD);

        //Fs
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        configStorageMSD.addListener(meterFs);

        //F
        DataSourceF meterF = new DataSourceF(configStorageMSD);
        configStorageMSD.addListener(meterF);

        //Percolation
        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorageMSD);
        atomFilterDeviation.setMinDistance(params.minDrFilter);
        atomFilterDeviation.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSD, atomFilterDeviation, params.log2StepS, params.log2StepE);
        configStorageMSD.addListener(meterPerc);

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

        //Pressure Tensor (G_inf)
        DataGroup dataPTensor = (DataGroup)pTensorAccumulator.getData();
        IData dataPTensorAvg = dataPTensor.getData(pTensorAccumulator.AVERAGE.index);
        IData dataPTensorSD = dataPTensor.getData(pTensorAccumulator.STANDARD_DEVIATION.index);
        double sd2PTensor = 0;
        for(int i=0; i<pIndex.length; i++){sd2PTensor += dataPTensorSD.getValue(pIndex[i])*dataPTensorSD.getValue(pIndex[i]);}
        sd2PTensor/=pIndex.length;


        String filenameVisc, filenameMSD, filenameD, filenameFs, filenameF, filenamePerc;

        if(sim.potentialChoice == SimGlass.PotentialChoice.HS){
            double phi;
            if(params.D == 2){
                phi = Math.PI/4*(params.nA+params.nB/(1.4*1.4))/volume;
            }else{
                phi = Math.PI/6*(params.nA+params.nB/(1.4*1.4*1.4))/volume;
            }
            System.out.println("rho: " + params.density + "  phi: " + phi+"\n");
            System.out.println("Z: " + pAvg/params.density/params.temperature +"  "+ pErr/params.density/params.temperature  +"  cor: "+pCorr);
            filenameVisc = String.format("viscRho%1.3f.out",rho);
            filenameMSD = String.format("msdRho%1.3f.out", rho);
            filenameD = String.format("dRho%1.3f.out", rho);
            filenameFs = String.format("fsRho%1.3f.out", rho);
            filenameF = String.format("fRho%1.3f.out", rho);
            filenamePerc = String.format("percRho%1.3f.out", rho);
        }else{
            DataGroup dataT = (DataGroup)tAccumulator.getData();
            IData dataTAvg = dataT.getData(tAccumulator.AVERAGE.index);
            IData dataTErr = dataT.getData(tAccumulator.ERROR.index);
            IData dataTCorr = dataT.getData(tAccumulator.BLOCK_CORRELATION.index);
            double tAvg  = dataTAvg.getValue(0);
            double tErr  = dataTErr.getValue(0);
            double tCorr = dataTCorr.getValue(0);
            System.out.println("rho: " + params.density+"\n");
            System.out.println("T: " + tAvg +"  "+ tErr +"  cor: "+tCorr);
            System.out.println("Z: " + pAvg/params.density/tAvg +"  "+ pErr/params.density/tAvg  +"  cor: "+pCorr);
            filenameVisc = String.format("viscPho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameMSD = String.format("msdPho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameD = String.format("dPho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameFs = String.format("fsPho%1.3fT%1.3f.out", rho, params.temperature);
            filenameF = String.format("fPho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePerc = String.format("percPho%1.3fT%1.3f.out",  rho, params.temperature);
            System.out.println("G: " + sim.box.getBoundary().volume()/tAvg*sd2PTensor+"\n");
        }
        System.out.println("P: " + pAvg +"  "+ pErr +"  cor: "+pCorr);

        //Viscosity
        FileWriter fileWriterVisc;
        try {
            fileWriterVisc = new FileWriter(filenameVisc, false);
            DataDoubleArray x = pTensorAccumVisc.getIndependentData(0);
            for (int i=0; i<pTensorAccumVisc.getData().getLength(); i++){
                double yi = pTensorAccumVisc.getData().getValue(i);
                double yiErr = pTensorAccumVisc.errData.getValue(i);
                double xi = x.getValue(i);
                if(!Double.isNaN(yi)){
                    fileWriterVisc.write(xi + " " + yi + " " + yiErr + "\n");
                }
            }
            fileWriterVisc.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }

        //MSD
        FileWriter fileWriterMSD, fileWriterD, fileWriterFs, fileWriterF, fileWriterPerc;
        try {
            fileWriterMSD = new FileWriter(filenameMSD,false);
            fileWriterD   = new FileWriter(filenameD,  false);
            fileWriterFs  = new FileWriter(filenameFs, false);
            fileWriterF   = new FileWriter(filenameF,  false);
            fileWriterPerc   = new FileWriter(filenamePerc,  false);
            DataDoubleArray x = meterMSD.getIndependentData(0);
            DataDoubleArray x2 = meterMSD.getIndependentData(0);
            for (int i=0; i<meterMSD.getData().getLength(); i++){
                double xi = x.getValue(i);
                double yi = meterMSD.getData().getValue(i);
                double yiErr = meterMSD.errData.getValue(i);
                double yiFs = meterFs.getData().getValue(i);
                double yiF  = meterF.getData().getValue(i);
                if( !Double.isNaN(yi) && !Double.isNaN(yiErr) && !Double.isNaN(yiFs)){
                    fileWriterMSD.write(xi + " " + yi + " " + yiErr+"\n");
                    fileWriterD.write(xi + " " + yi/6/xi + " " + yiErr/6/xi + "\n");
                    fileWriterFs.write(xi + " " + yiFs + "\n");
                    fileWriterF.write(xi + " " + yiF + "\n");
                    if(i >= params.log2StepS && i <= params.log2StepE){
                        double yiPerc  = meterPerc.getData().getValue(i);
                        fileWriterPerc.write(xi + " " + yiPerc + "\n");
                    }
                }
            }
            fileWriterMSD.close();
            fileWriterD.close();
            fileWriterFs.close();
            fileWriterF.close();
            fileWriterPerc.close();
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