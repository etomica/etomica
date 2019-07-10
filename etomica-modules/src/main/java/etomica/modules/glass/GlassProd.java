/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.*;
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
            params.nA = 125;
            params.nB = 125;
            params.density = 1.25;
            params.D = 3;
            params.temperature = 1;
            params.numStepsEq = 1000;
            params.numSteps = 10000;
            params.minDrFilter = 0.4;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);
        System.out.println(params.D +"D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        double volume = sim.box.getBoundary().volume();
        int numAtoms = params.nA + params.nB;
        double rho= numAtoms/volume;
        System.out.println("T = " + params.temperature);
        System.out.println( params.numSteps + " MD steps after " + params.numStepsEq + " equilibaration steps");

        //Equilibration
        double temperature0 = params.temperatureMelt > params.temperature ? params.temperatureMelt : params.temperature;
        if (temperature0 > params.temperature) System.out.println("Equilibrating at T=" + temperature0);
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(temperature0);
        sim.activityIntegrate.setMaxSteps(params.numStepsEq);
        sim.getController().actionPerformed();
        sim.getController().reset();

        if (temperature0 > params.temperature) {
            System.out.println("Equilibrating at T=" + params.temperature);
            sim.integrator.setTemperature(params.temperature);
            sim.getController().actionPerformed();
            sim.getController().reset();
        }

        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);
        sim.activityIntegrate.setMaxSteps(params.numSteps);

        long blocksize = params.numSteps / 100;

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

        //Viscosity
        DataFork pTensorFork = new DataFork();
        AccumulatorPTensor pTensorAccumVisc = new AccumulatorPTensor(sim.integrator, sim.integrator.getTimeStep());
        pTensorAccumVisc.setEnabled(true);
        pTensorFork.addDataSink(pTensorAccumVisc);
        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork);


        AccumulatorAverageFixed pTensorAccumulator = new AccumulatorAverageFixed(1);
        pTensorFork.addDataSink(pTensorAccumulator);

        sim.integrator.getEventManager().addListener(pTensorAccumViscPump);

        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        pFork.addDataSink(pAccumulator);

        AccumulatorAverageFixed tAccumulator = null;
        AccumulatorAverageFixed accPE = null;

        if(sim.potentialChoice != SimGlass.PotentialChoice.HS){
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotentialMaster(), sim.box);
            long bs = blocksize / 5;
            if (bs == 0) bs = 1;
            accPE = new AccumulatorAverageFixed(bs);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 5);
            sim.integrator.getEventManager().addListener(pumpPE);

            MeterTemperature tMeter = new MeterTemperature(sim, sim.box, params.D);
            tAccumulator = new AccumulatorAverageFixed(blocksize / 5);
            DataPumpListener tPump = new DataPumpListener(tMeter, tAccumulator, 5);
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
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSD, atomFilterDeviation, 5, 30);

        //Immobile fraction
        DataSourcePercolation.ImmFractionSource meterImmFraction = meterPerc.makeImmFractionSource();

        configStorageMSD.addListener(meterPerc);

        //Strings
        DataSourceStrings meterL = new DataSourceStrings(configStorageMSD, 3, 30);
        configStorageMSD.addListener(meterL);


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

        String filenameVisc, filenameMSD, filenameD, filenameFs, filenameF, filenamePerc, filenameImmFrac, filenameL;

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
            filenameL = String.format("lRho%1.3f.out", rho);
            filenameImmFrac = String.format("immFracRho%1.3f.out", rho);
        }else{
            // Energy
            DataGroup dataU = (DataGroup) accPE.getData();
            IData dataUAvg = dataU.getData(accPE.AVERAGE.index);
            IData dataUErr = dataU.getData(accPE.ERROR.index);
            IData dataUCorr = dataU.getData(accPE.BLOCK_CORRELATION.index);
            double uAvg = dataUAvg.getValue(0);
            double uErr = dataUErr.getValue(0);
            double uCorr = dataUCorr.getValue(0);

            //Pressure Tensor (G_inf)
            DataGroup dataPTensor = (DataGroup) pTensorAccumulator.getData();
            IData dataPTensorSD = dataPTensor.getData(pTensorAccumulator.STANDARD_DEVIATION.index);
            double sd2PTensor = 0;
            for (int i = 0; i < pIndex.length; i++) {
                sd2PTensor += dataPTensorSD.getValue(pIndex[i]) * dataPTensorSD.getValue(pIndex[i]);
            }
            sd2PTensor /= pIndex.length;

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
            System.out.println("U: " + uAvg / numAtoms + "  " + uErr / numAtoms + "  cor: " + uCorr);
            filenameVisc = String.format("viscRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameMSD = String.format("msdRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameD = String.format("dRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameFs = String.format("fsRho%1.3fT%1.3f.out", rho, params.temperature);
            filenameF = String.format("fRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePerc = String.format("percRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameL = String.format("lRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameImmFrac = String.format("immFracRho%1.3fT%1.3f.out",  rho, params.temperature);
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
        FileWriter fileWriterMSD, fileWriterD, fileWriterFs, fileWriterF, fileWriterPerc, fileWriterImmFrac, fileWriterL;
        try {
            fileWriterMSD = new FileWriter(filenameMSD,false);
            fileWriterD   = new FileWriter(filenameD,  false);
            fileWriterFs  = new FileWriter(filenameFs, false);
            fileWriterF   = new FileWriter(filenameF,  false);
            fileWriterPerc   = new FileWriter(filenamePerc,  false);
            fileWriterImmFrac   = new FileWriter(filenameImmFrac ,  false);
            fileWriterL   = new FileWriter(filenameL ,  false);
            DataDoubleArray x = meterMSD.getIndependentData(0);
            for (int i=0; i<meterMSD.getData().getLength(); i++){
                double xi = x.getValue(i);
                double yi = meterMSD.getData().getValue(i);
                double yiErr = meterMSD.errData.getValue(i);
                double yiFs = meterFs.getData().getValue(i);
                double yiF  = meterF.getData().getValue(i);
                double yiL  = meterL.getData().getValue(i);
                if( !Double.isNaN(yi) && !Double.isNaN(yiErr) && !Double.isNaN(yiFs)){
                    fileWriterMSD.write(xi + " " + yi + " " + yiErr+"\n");
                    fileWriterD.write(xi + " " + yi/6/xi + " " + yiErr/6/xi + "\n");
                    fileWriterFs.write(xi + " " + yiFs + "\n");
                    fileWriterF.write(xi + " " + yiF + "\n");
                    double yiPerc  = meterPerc.getData().getValue(i);
                    double yiImmFrac  = meterImmFraction.getData().getValue(i);
                    if(!Double.isNaN(yiPerc)){
                        fileWriterPerc.write(xi + " " + yiPerc + "\n");
                    }
                    if(!Double.isNaN(yiImmFrac)){
                        fileWriterImmFrac.write(xi + " " + yiImmFrac + "\n");
                    }
                    if(!Double.isNaN(yiL)){
                        fileWriterL.write(xi + " " + yiL + "\n");
                    }
                }
            }
            fileWriterMSD.close();
            fileWriterD.close();
            fileWriterFs.close();
            fileWriterF.close();
            fileWriterPerc.close();
            fileWriterImmFrac.close();
            fileWriterL.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        double time1 = System.currentTimeMillis();
        System.out.println("\ntime: " + (time1-time0)/1000/3600 + " hrs");
    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;
        public double minDrFilter = 0.4;
        public double temperatureMelt = 0;
    }
}