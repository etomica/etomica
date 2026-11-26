/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.*;
import etomica.data.meter.*;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;

/**
 * This class normally consists of 3 stages:
 * 1. NVT equilibration (numStepsEqIsothermal)
 * 2. NVT collection of total energy (numStepsIsothermal)
 * 3. NVE production run at total energy measured in stage 2 (numSteps)
 *
 * For a soft potential that has trouble melting from the initial configuration, 4 stages can be used:
 * 1. NVT equilibration at a temperature high enough to melt (numStepsEqIsothermal)
 * 2. NVT equilibration at temperature of interest (numStepsEqIsothermal)
 * 3. NVT collection of total energy at temperature of interest (numStepsIsothermal)
 * 4. NVE production run at total energy measured in stage 3 (numSteps)
 *
 * After a simulation finishes, the simulation can be continued from where it ended.  This allows a very long run (which
 * may be necessary to collect long-time dynamics) can be split up across many jobs.  In each case, the run should be
 * specified as consisting of one or more of the 3 stages described above.  A stage which should not be run should have
 * its steps set to 0.  Save/restore does not understand the 4-stage approach; if needed, just do several equilibration
 * runs with different temperatures.  Whatever type of simulation, the job needs to finish (the class does not handle
 * resuming after an interruption).  Each job can continue the last stage from the previous job, or can start the next
 * stage.
 */

public class TransportCoefficients {

    public static void main(String[] args) {
        long time0eq = System.nanoTime();
        double startTime = System.nanoTime()/1e9;
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 100;
            params.nB = 100;
            params.density =  0.5; // 2D params.density = 0.509733; //3D  = 0.69099;
            params.D = 3;
            params.temperature = 1.0;
            params.numStepsEqIsothermal = 10000;
            params.numStepsIsothermal = 0;
            params.numSteps = 100000;
            params.rcLJ = 2.5;
            params.tStep = 0.01;
            params.numStepsIsothermal = 0;
        }
        if (params.nB==0) params.doSwap = false;

        double rho = params.density;
        if (params.eta>0 && params.potential == SimGlass.PotentialChoice.HS) {
            rho = SimGlass.densityForEta(params.eta, params.nA/(double)(params.nA+params.nB), params.sigmaB);
            params.density = rho;
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep, params.rcLJ, params.sigmaB);
        System.out.println(params.D + "D " + sim.potentialChoice);
        System.out.println("L: " + sim.box.getBoundary().getBoxSize().getX(0));
        System.out.println("nA:nB " + params.nA + ":" + params.nB);
        if (params.nB > 0 ) {
            System.out.println("sigmaA: " + 1.0 + " sigmaB: " + sim.sigmaB);
        }

        int numAtoms = params.nA + params.nB;
        System.out.println("T = " + params.temperature);
        System.out.println("doSwap: " + params.doSwap);
        System.out.println(params.numSteps + " production MD steps after " + params.numStepsEqIsothermal + " NVT equilibration steps, "+params.numStepsIsothermal+" NVT E check using dt = " + sim.integrator.getTimeStep());
        double volume = sim.box.getBoundary().volume();
        if (params.potential == SimGlass.PotentialChoice.HS) {
            double phi;
            double sb = params.sigmaB == 0 ? 1.0/1.4 : params.sigmaB;
            if (params.nB == 0) sb = 1.0;
            if (params.D == 2) {
                phi = Math.PI / 4 * (params.nA + params.nB * sb * sb) / volume;
            } else {
                phi = Math.PI / 6 * (params.nA + params.nB * sb * sb * sb) / volume;
            }
            System.out.println("rho: " + params.density + "  phi: " + phi + "\n");
        } else {
            System.out.println("rho: " + params.density + "\n");
        }
        sim.initConfig();

        long blocksize = params.numSteps / 100;

        //Equilibration
        double temperature0 = Math.max(params.temperatureMelt, params.temperature);
        if (temperature0 > params.temperature) System.out.println("Equilibrating at T=" + temperature0);
        sim.integrator.setIsothermal(true);
        if (params.doSwap) sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(temperature0);
        // Equilibrate isothermally
        double maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsEqIsothermal, false, maxTime));
        if (System.nanoTime()/1e9 - startTime >= params.maxWalltime) {
            // we've reached our time limit.  disable the rest of the sim.
            params.numSteps = params.numStepsIsothermal = 0;
            // and set our eq steps to what we were able to actually run
            System.out.println(sim.integrator.getStepCount()+" steps of isothermal equilibration finished before time limit");
        }

        AccumulatorAverageFixed accE = new AccumulatorAverageFixed(1);
        MeterEnergyFromIntegrator meterE = new MeterEnergyFromIntegrator(sim.integrator);
        DataPumpListener pumpE = new DataPumpListener(meterE, accE, 10);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            sim.integrator.getEventManager().addListener(pumpE);
        }


        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            if (temperature0 == params.temperature) {
                // Now collect average total energy (stage 2).
                maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
                sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsIsothermal, false, maxTime));
                if (System.nanoTime()/1e9 - startTime >= params.maxWalltime) {
                    // we've reached our time limit.  disable the rest of the sim.
                    params.numSteps = 0;
                    // and set our nvt steps to what we were able to actually run
                    System.out.println(sim.integrator.getStepCount()+" steps of isothermal finished before time limit");
                }
            }
            else {
                // Don't be trying to save/restore into here.
                // Just do a run at a high T, and then another at the T you want.

                // first stage was at higher temperature to ensure crystal melted
                // now we need a second stage at the temperature we actually want
                System.out.println("Equilibrating at T=" + params.temperature);
                sim.integrator.setTemperature(params.temperature);
                sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsEqIsothermal));
                accE.reset();

                // Assume isothermal equilibration.  Now collect average total energy.
                sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsIsothermal));
            }
        }

        long time2eq = System.nanoTime();
        if (params.numStepsEqIsothermal + params.numStepsIsothermal > 0 && params.doSwap) {
            // only for this sim, swapMove not saved/restored (has no state, only tracker has state)
            System.out.println("swap acceptance: "+sim.swapMove.getTracker().acceptanceProbability());
        }
        System.out.printf("equilibration time: %3.2f min\n", (time2eq-time0eq)/1e9/60.);

        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            // Try to use the average total energy from the 2nd half of equilibration
            // to set the current energy now (via the temperature) so that the average
            // temperature during production will be approximately equal to the set
            // temperature.
            double avgE = accE.getData(accE.AVERAGE).getValue(0);
            System.out.println("average energy during second half of eq: " + avgE / numAtoms);
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            double pe = meterPE.getDataAsScalar();
//            double ke = new MeterKineticEnergy(sim.box).getDataAsScalar();
//            System.out.println("kinetic energy at end of eq: "+ke/numAtoms);
//            System.out.println("potential energy at end of eq: "+pe/numAtoms);
            double newKE = avgE - pe;
            double nowTemp = newKE * 2.0 / params.D / numAtoms;
            double oldTemp = new MeterTemperature(sim.box, params.D).getDataAsScalar();
            System.out.println("setting temp " + oldTemp + " => " + nowTemp); // + " (ke "+ke+" => " + newKE / numAtoms + ")");
            if (params.doSwap) sim.integrator.setIntegratorMC(null, 0);
            sim.integrator.setIsothermal(false);
            sim.integrator.setTemperature(nowTemp);
            accE.reset();
//            double newNewKE = new MeterKineticEnergy(sim.box).getDataAsScalar();
//            System.out.println("actual KE => "+newNewKE/numAtoms);
        }


        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);

        //RDF-A
        final MeterRDF meterRdfA = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListenerA = new IntegratorListenerAction(meterRdfA);
        sim.integrator.getEventManager().addListener(rdfMeterListenerA);
        rdfMeterListenerA.setInterval(100);
        meterRdfA.setAtomType(sim.speciesA.getLeafType());
        meterRdfA.getXDataSource().setXMax(1.5);
        meterRdfA.setBox(sim.box);


        //RDF-B
        final MeterRDF meterRdfB = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListenerB = new IntegratorListenerAction(meterRdfB);
        sim.integrator.getEventManager().addListener(rdfMeterListenerB);
        rdfMeterListenerB.setInterval(100);
        meterRdfB.setAtomType(sim.speciesB.getLeafType());
        meterRdfB.getXDataSource().setXMax(sim.sigmaB * 1.5);
        meterRdfB.setBox(sim.box);


        //RDF-AB
        final MeterRDF meterRdf = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListener = new IntegratorListenerAction(meterRdf);
        sim.integrator.getEventManager().addListener(rdfMeterListener);
        rdfMeterListener.setInterval(100);
        double sigmaAB = 1.0;
        meterRdf.setAtomType(sim.speciesA.getLeafType());
        if (params.nB > 0) {
            meterRdf.setAtomTypes(sim.speciesA.getLeafType(), sim.speciesB.getLeafType());
            sigmaAB = (1.0 + sim.sigmaB) / 2;
        }
        meterRdf.getXDataSource().setXMax(sigmaAB * 1.5);
        meterRdf.setBox(sim.box);


        //P
        IDataSource pTensorMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pTensorMeter = new MeterPressureTensor(sim.integrator.getPotentialCompute(), sim.box);
        } else {
            pTensorMeter = new MeterPressureHardTensor((IntegratorHard) sim.integrator);
            ((MeterPressureHardTensor)pTensorMeter).setDoNonEquilibrium(true);
        }

        //Viscosity
        DataFork pTensorFork = new DataFork();
        int dn = 1;
        AccumulatorPTensor pTensorAccumVisc = new AccumulatorPTensor(sim.box, dn * sim.integrator.getTimeStep());
        pTensorAccumVisc.setEnabled(true);
        pTensorFork.addDataSink(pTensorAccumVisc);

        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork, dn);


        sim.integrator.getEventManager().addListener(pTensorAccumViscPump);

        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        DataFork pFork = new DataFork();
        tracer.setDataSink(pFork);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        pFork.setDataSink(pAccumulator);

        AccumulatorAverageFixed tAccumulator = null;
        AccumulatorAverageFixed accPE = null;

        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            long bs = blocksize / 5;
            if (bs == 0) bs = 1;
            accPE = new AccumulatorAverageFixed(bs);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 5);
            sim.integrator.getEventManager().addListener(pumpPE);

            MeterTemperature tMeter = new MeterTemperature(sim.getSpeciesManager(), sim.box, params.D);
            tAccumulator = new AccumulatorAverageFixed(blocksize / 5);
            DataPumpListener tPump = new DataPumpListener(tMeter, tAccumulator, 5);
            tAccumulator.addDataSink(pTensorAccumVisc.makeTemperatureSink(), new AccumulatorAverage.StatType[]{tAccumulator.AVERAGE});
            sim.integrator.getEventManager().addListener(tPump);
        }


         //MSD
        ConfigurationStorage configStorageMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
        configStorageMSD.setEnabled(true);
        ConfigurationStorage configStorageMSD3 = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD, 60, 3);
        configStorageMSD3.setEnabled(true);
        DataSourceMSD meterMSD = new DataSourceMSD(configStorageMSD);
        configStorageMSD.addListener(meterMSD);

        DataSourceMSD meterMSDA = new DataSourceMSD(configStorageMSD, sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterMSDA);
        DataSourceMSD meterMSDB = new DataSourceMSD(configStorageMSD, sim.speciesB.getLeafType());
        if (params.nB > 0) configStorageMSD.addListener(meterMSDB);

        AtomStressSource stressSource;
        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            AtomHardStressCollector ahsc = new AtomHardStressCollector(sim.integrator);
            ((IntegratorHard) sim.integrator).addCollisionListener(ahsc);
            stressSource = ahsc;
        } else {
            stressSource = new PotentialCallbackAtomStress(sim.box);
        }
        int[][] normalComps = new int[sim.getSpace().getD()][2];
        for (int i = 0; i < normalComps.length; i++) {
            normalComps[i][0] = normalComps[i][1] = i;
        }

        sim.integrator.getEventManager().addListener(configStorageMSD3);
        sim.integrator.getEventManager().addListener(configStorageMSD);


        //Production ...
        long time0 = System.nanoTime();
        maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps, false, maxTime));
        params.numSteps = sim.integrator.getStepCount();
        System.out.println(params.numSteps+" steps of production finished before time limit");

        //Pressure
        DataGroup dataP = (DataGroup) pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg = dataPAvg.getValue(0);
        double pErr = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);

        String filenameVisc, filenameMSD, filenameMSDA, filenameMSDB, filenameRdfA, filenameRdfB, filenameRdf;

        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            System.out.println("Z: " + pAvg / params.density / params.temperature + "  " + pErr / params.density / params.temperature + "  cor: " + pCorr);
            filenameVisc = String.format("viscRho%1.3f.dat", rho);
            filenameMSD = String.format("msdRho%1.3f.dat", rho);
            filenameMSDA = String.format("msdARho%1.3f.dat", rho);
            filenameMSDB = String.format("msdBRho%1.3f.dat", rho);
            filenameRdfA = String.format("rdfARho%1.3f.dat", rho);
            filenameRdfB = String.format("rdfBRho%1.3f.dat", rho);
            filenameRdf = String.format("rdfRho%1.3f.dat", rho);
        } else {
            // Energy
            DataGroup dataU = (DataGroup) accPE.getData();
            IData dataUAvg = dataU.getData(accPE.AVERAGE.index);
            IData dataUErr = dataU.getData(accPE.ERROR.index);
            IData dataUCorr = dataU.getData(accPE.BLOCK_CORRELATION.index);
            double uAvg = dataUAvg.getValue(0);
            double uErr = dataUErr.getValue(0);
            double uCorr = dataUCorr.getValue(0);

            DataGroup dataT = (DataGroup) tAccumulator.getData();
            IData dataTAvg = dataT.getData(tAccumulator.AVERAGE.index);
            IData dataTErr = dataT.getData(tAccumulator.ERROR.index);
            IData dataTCorr = dataT.getData(tAccumulator.BLOCK_CORRELATION.index);
            double tAvg = dataTAvg.getValue(0);
            double tErr = dataTErr.getValue(0);
            double tCorr = dataTCorr.getValue(0);

            IData dataE = accE.getData();
            double eAvg = dataE.getValue(accE.AVERAGE.index);
            double eStdev = dataE.getValue(accE.STANDARD_DEVIATION.index);
            double eCor = dataE.getValue(accE.BLOCK_CORRELATION.index);

            System.out.println("T: " + tAvg + "  " + tErr + "  cor: " + tCorr);
            System.out.println("Z: " + pAvg / params.density / tAvg + "  " + pErr / params.density / tAvg + "  cor: " + pCorr);
            System.out.println("U: " + uAvg / numAtoms + "  " + uErr / numAtoms + "  cor: " + uCorr);
            System.out.println("E: " + eAvg / numAtoms + "  " + eStdev / numAtoms + "  cor: " + eCor);
            filenameVisc = String.format("viscRho%1.3fT%1.3f.dat", rho, params.temperature);

            filenameMSD = String.format("msdRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameMSDA = String.format("msdARho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameMSDB = String.format("msdBRho%1.3fT%1.3f.dat", rho, params.temperature);

            filenameRdfA = String.format("rdfARho%1.3f.dat", rho);
            filenameRdfB = String.format("rdfBRho%1.3f.dat", rho);
            filenameRdf = String.format("rdfRho%1.3f.dat", rho);

        }
        System.out.println("P: " + pAvg + "  " + pErr + "  cor: " + pCorr);

        try {
            TransportCoefficients.writeDataToFile(pTensorAccumVisc, pTensorAccumVisc.errData, filenameVisc);
            TransportCoefficients.writeDataToFile(meterMSD, meterMSD.getError(), filenameMSD);
            TransportCoefficients.writeDataToFile(meterRdf, filenameRdf);
            if (params.nB > 0) {
                TransportCoefficients.writeDataToFile(meterMSDA, meterMSDA.getError(), filenameMSDA);
                TransportCoefficients.writeDataToFile(meterMSDB, meterMSDB.getError(), filenameMSDB);
                TransportCoefficients.writeDataToFile(meterRdfA, filenameRdfA);
                TransportCoefficients.writeDataToFile(meterRdfB, filenameRdfB);
            }
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        long time1 = System.nanoTime();
        double sim_time = (time1 - time0) / 1e9 / 60.0;
        System.out.printf("\ntime: %3.2f min\n", sim_time);
    }

    public static class SimParams extends SimGlass.GlassParams {
        public long numStepsEqIsothermal = 10000;
        public long numSteps = 1000000;
        public double temperatureMelt = 0;
        public double maxWalltime = Double.POSITIVE_INFINITY;
        public long numStepsIsothermal = 0;
    }

    public static void writeDataToFile(IDataSource meter, String filename) throws IOException {
        writeDataToFile(meter, null, filename);
    }

    public static void writeDataToFile(IDataSource meter, IData data2, String filename) throws IOException {
        writeDataToFile(meter, data2, null, filename);
    }

    public static void writeDataToFile(IDataSource meter, IData data2, IData data3, String filename) throws IOException {
        IData data;
        IData xData;
        if (meter instanceof AccumulatorAverage) {
            AccumulatorAverage acc = (AccumulatorAverage) meter;
            data = acc.getData(acc.AVERAGE);
            xData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) acc.getDataInfo()).getSubDataInfo(acc.AVERAGE.index)).getXDataSource().getIndependentData(0);
        } else {
            data = meter.getData();
            xData = ((DataFunction.DataInfoFunction) meter.getDataInfo()).getXDataSource().getIndependentData(0);
        }
        boolean allNaN = true;
        for (int i = 0; i < xData.getLength(); i++) {
            if (!Double.isNaN(data.getValue(i))) allNaN = false;
        }
        if (allNaN) return;
        FileWriter fw = new FileWriter(filename);
        for (int i = 0; i < xData.getLength(); i++) {
            double y = data.getValue(i);
            if (Double.isNaN(y)) continue;
            if (data2 == null) {
                fw.write(xData.getValue(i) + " " + y + "\n");
            }
            else if (data3 == null) {
                fw.write(xData.getValue(i) + " " + y + " " + data2.getValue(i) + "\n");
            }
            else {
                fw.write(xData.getValue(i) + " " + y + " " + data2.getValue(i) + " " +data3.getValue(i) + "\n");
            }
        }
        fw.close();
    }

}