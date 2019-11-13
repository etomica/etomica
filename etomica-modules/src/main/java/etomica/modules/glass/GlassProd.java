/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
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
            params.nA = 50;
            params.nB = 50;
            params.density = 1.2; // 2D params.density = 0.509733; //3D  = 0.69099;
            params.D = 3;
            params.temperature = 1.0;
            params.numStepsEq = 100000;
            params.numSteps =   1000000;
            params.minDrFilter = 0.4;
            params.qx = 7.0;
            params.log2StepMin = 10;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);
        System.out.println(params.D +"D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        double volume = sim.box.getBoundary().volume();
        int numAtoms = params.nA + params.nB;
        double rho= numAtoms/volume;
        System.out.println("T = " + params.temperature);
        System.out.println( params.numSteps + " MD steps after " + params.numStepsEq + " equilibaration steps , using dt = " + sim.integrator.getTimeStep());

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
        int dn = 1;
        AccumulatorPTensor pTensorAccumVisc = new AccumulatorPTensor(sim.integrator, dn*sim.integrator.getTimeStep());
        pTensorAccumVisc.setEnabled(true);
        pTensorFork.addDataSink(pTensorAccumVisc);
        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork, dn);

        DataProcessor dpSquared = new DataProcessor() {
            DataDouble data = new DataDouble();

            @Override
            protected IData processData(IData inputData) {
                data.x = 0;
                int n = 0;
                for (int i = 0; i < params.D; i++) {
                    for (int j = i + 1; j < params.D; j++) {
                        double c = ((DataTensor) inputData).x.component(i, j);
                        data.x += c * c;
                        n++;
                    }
                }
                data.x /= n;
                return data;
            }

            @Override
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                dataInfo = new DataDouble.DataInfoDouble("G", Null.DIMENSION);
                return dataInfo;
            }
        };
        pTensorFork.addDataSink(dpSquared);

        AccumulatorAverageFixed gTensorAccumulator = new AccumulatorAverageFixed(blocksize);
        dpSquared.setDataSink(gTensorAccumulator);

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


        //shear stress AC
        AccumulatorAutocorrelationShearStress dpxyAutocor = new AccumulatorAutocorrelationShearStress(256, sim.integrator.getTimeStep());
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            pTensorFork.addDataSink(dpxyAutocor);
        }
        dpxyAutocor.setPushInterval(16384);


        //MSD
        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSD.setEnabled(true);
        DataSourceMSD meterMSD = new DataSourceMSD(configStorageMSD);
        configStorageMSD.addListener(meterMSD);

        //VAC
        DataSourceVAC meterVAC = new DataSourceVAC(configStorageMSD);
        configStorageMSD.addListener(meterVAC);

        //Fs
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        Vector q = sim.getSpace().makeVector();
        q.setX(0, params.qx);
        meterFs.setQ(q);
        configStorageMSD.addListener(meterFs);

        //F
//        DataSourceF meterF = new DataSourceF(configStorageMSD);
//        configStorageMSD.addListener(meterF);

        //Percolation
        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorageMSD);
        atomFilterDeviation.setMinDistance(params.minDrFilter);
        atomFilterDeviation.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSD, atomFilterDeviation, params.log2StepMin, 30);

        //Immobile fraction
        DataSourcePercolation.ImmFractionSource meterImmFraction = meterPerc.makeImmFractionSource();

        configStorageMSD.addListener(meterPerc);

        DataSourcePercolation0 meterPerc0 = new DataSourcePercolation0(sim.box, sim.getRandom());
        meterPerc0.setImmFracs(new double[]{0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                0.70, 0.75, 0.8, 0.85, 0.9, 0.95, 1});
        AccumulatorAverageFixed accPerc0 = new AccumulatorAverageFixed(10);
        DataPumpListener pumpPerc0 = new DataPumpListener(meterPerc0, accPerc0, 1000);
        sim.integrator.getEventManager().addListener(pumpPerc0);

        //Strings
        DataSourceStrings meterL = new DataSourceStrings(configStorageMSD, 3, 30);
        configStorageMSD.addListener(meterL);

        //Alpha2
        DataSourceAlpha2 meterAlpha2 = new DataSourceAlpha2(configStorageMSD);
        configStorageMSD.addListener(meterAlpha2);


        //S(q)

        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);  // just average, no uncertainty
        accSFac.setPushInterval(1);
        MeterStructureFactor meterSFac = new MeterStructureFactor(sim.getSpace(), sim.box, 15);
        DataPumpListener pumpSFac = new DataPumpListener(meterSFac, accSFac, 1000);
        sim.integrator.getEventManager().addListener(pumpSFac);
        double vB = sim.getSpace().powerD(sim.sigmaB);
        meterSFac.setAtomTypeFactor(sim.speciesB.getAtomType(0), vB);


        sim.integrator.getEventManager().addListener(configStorageMSD);

        //Run
        double time0 = System.currentTimeMillis();
        sim.getController().actionPerformed();

        //Pressure
        DataGroup dataP = (DataGroup)pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg  = dataPAvg.getValue(0);
        double pErr  = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);


        String filenameVisc, filenameMSD, filenameD, filenameFs, filenameF, filenamePerc, filenamePerc0, filenameImmFrac, filenameL, filenameAlpha2, filenameSq, filenameVAC;

        String filenamePxyAC = "";
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
            filenameVAC = String.format("vacRho%1.3f.out", rho);
            filenameFs = String.format("fsRho%1.3fQ%1.2f.out", rho, params.qx);
//            filenameF = String.format("fRho%1.3f.out", rho);
            filenamePerc = String.format("percRho%1.3f.out", rho);
            filenamePerc0 = String.format("perc0Rho%1.3f.out", rho);
            filenameL = String.format("lRho%1.3f.out", rho);
            filenameImmFrac = String.format("immFracRho%1.3f.out", rho);
            filenameAlpha2 = String.format("alpha2Rho%1.3f.out", rho);
            filenameSq = String.format("sqRho%1.3f.out",rho);
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
            IData dataG = gTensorAccumulator.getData();
            double avgG = dataG.getValue(gTensorAccumulator.AVERAGE.index);
            double errG = dataG.getValue(gTensorAccumulator.ERROR.index);
            double corG = dataG.getValue(gTensorAccumulator.BLOCK_CORRELATION.index);

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
            filenameVAC = String.format("vacRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameFs = String.format("fsRho%1.3fT%1.3fQ%1.2f.out", rho, params.temperature, params.qx);
//            filenameF = String.format("fRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePerc = String.format("percRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePerc0 = String.format("perc0Rho%1.3fT%1.3f.out", rho, params.temperature);
            filenameL = String.format("lRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameImmFrac = String.format("immFracRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameAlpha2 = String.format("alpha2Rho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePxyAC = String.format("acPxyRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameSq = String.format("sqRho%1.3fT%1.3f.out",  rho, params.temperature);
            double V = sim.box.getBoundary().volume();
            System.out.println("G: " + V * avgG / tAvg + " " + V * errG / tAvg + " cor: " + corG + "\n");
        }
        System.out.println("P: " + pAvg +"  "+ pErr +"  cor: "+pCorr);

        // shear stress AC
        if(sim.potentialChoice != SimGlass.PotentialChoice.HS){
            try {
                FileWriter fileWriterPxyAC = new FileWriter(filenamePxyAC, false);
                DataDoubleArray x = dpxyAutocor.getIndependentData(0);
                for (int i=0; i<dpxyAutocor.getData().getLength(); i++){
                    double yi = dpxyAutocor.getData().getValue(i);
                    double xi = x.getValue(i);
                    if(!Double.isNaN(yi)){
                        fileWriterPxyAC.write(xi + " " + yi + "\n");
                    }
                }
                fileWriterPxyAC.close();
            } catch (IOException e) {
                System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
            }
        }


        //S(q)
        DataGroup dataSF = (DataGroup)accSFac.getData();
        IData dataSFAvg = dataSF.getData(accSFac.AVERAGE.index);
        int nSF  = dataSFAvg.getLength();
        IData xData = meterSFac.getIndependentData(0);
        try {
            FileWriter fileWriterSq = new FileWriter(filenameSq, false);
            for(int i=0;i<nSF; i++){
                fileWriterSq.write(xData.getValue(i) + " "+ dataSFAvg.getValue(i)+"\n");
            }
            fileWriterSq.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }




        //Viscosity
        try {
            FileWriter fileWriterVisc = new FileWriter(filenameVisc, false);
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
        try {
            FileWriter fileWriterMSD = new FileWriter(filenameMSD, false);
            FileWriter fileWriterD = new FileWriter(filenameD, false);
            FileWriter fileWriterVAC = new FileWriter(filenameVAC, false);
            FileWriter fileWriterFs = new FileWriter(filenameFs, false);
//           FileWriter fileWriterF   = new FileWriter(filenameF,  false);
            FileWriter fileWriterPerc = new FileWriter(filenamePerc, false);
            FileWriter fileWriterImmFrac = new FileWriter(filenameImmFrac, false);
            FileWriter fileWriterL = new FileWriter(filenameL, false);
            FileWriter fileWriterAlpha2 = new FileWriter(filenameAlpha2, false);
            DataDoubleArray x = meterMSD.getIndependentData(0);
            for (int i=0; i<meterMSD.getData().getLength(); i++){
                double xi = x.getValue(i);
                double yi = meterMSD.getData().getValue(i);
                double yiErr = meterMSD.errData.getValue(i);
                double yiFs = meterFs.getData().getValue(i);
//                double yiF  = meterF.getData().getValue(i);
                double yiVAC = meterVAC.getData().getValue(i);
                double yiVACErr = meterVAC.errData.getValue(i);
                double yiL  = meterL.getData().getValue(i);
                if( !Double.isNaN(yi) && !Double.isNaN(yiErr) && !Double.isNaN(yiFs)){
                    fileWriterMSD.write(xi + " " + yi + " " + yiErr+"\n");
                    fileWriterD.write(xi + " " + yi/2/params.D/xi + " " + yiErr/2/params.D/xi + "\n");
                    fileWriterVAC.write(xi + " " + yiVAC + " " + yiVACErr +"\n");
                    fileWriterFs.write(xi + " " + yiFs + "\n");
//                    fileWriterF.write(xi + " " + yiF + "\n");
                    double yiPerc  = meterPerc.getData().getValue(i);
                    double yiImmFrac  = meterImmFraction.getData().getValue(i);
                    double yiAlpha2  = meterAlpha2.getData().getValue(i);

                    if(!Double.isNaN(yiPerc)){
                        fileWriterPerc.write(xi + " " + yiPerc + "\n");
                    }
                    if(!Double.isNaN(yiImmFrac)){
                        fileWriterImmFrac.write(xi + " " + yiImmFrac + "\n");
                    }
                    if(!Double.isNaN(yiL)){
                        fileWriterL.write(xi + " " + yiL + "\n");
                    }
                    if(!Double.isNaN(yiAlpha2)){
                        fileWriterAlpha2.write(xi + " " + yiAlpha2 + "\n");
                    }

                }
            }
            fileWriterMSD.close();
            fileWriterVAC.close();
            fileWriterD.close();
            fileWriterFs.close();
//            fileWriterF.close();
            fileWriterPerc.close();
            fileWriterImmFrac.close();
            fileWriterL.close();
            fileWriterAlpha2.close();
            FileWriter fileWriterPerc0 = new FileWriter(filenamePerc0, false);
            x = meterPerc0.getIndependentData(0);
            IData perc0Avg = accPerc0.getData(accPerc0.AVERAGE);
            for (int i = 0; i < x.getLength(); i++) {
                double xi = x.getValue(i);
                double yi = perc0Avg.getValue(i);
                fileWriterPerc0.write(xi + " " + yi + "\n");
            }
            fileWriterPerc0.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        double time1 = System.currentTimeMillis();
        double sim_time = (time1-time0)/1000/3600;
        if(sim_time > 1.0){
            System.out.println("\ntime: " + sim_time + " hrs");
        }else{
            System.out.println("\ntime: " + (sim_time*60) + " mins");
        }


    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;
        public double minDrFilter = 0.4;
        public int log2StepMin = 5;
        public double temperatureMelt = 0;
        public double qx = 7.0;
    }
}
