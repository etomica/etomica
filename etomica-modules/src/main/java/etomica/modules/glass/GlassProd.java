/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.ParseArgs;
import etomica.util.Statefull;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

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

public class GlassProd {

    public static void saveObjects(StepsState steps, List<Statefull> objects) {
        try {
            FileWriter fw = new FileWriter("glass.steps");
            steps.saveState(fw);
            fw.close();
            fw = new FileWriter("glass.state");
            for (Statefull s : objects) {
                if (s!=null) s.saveState(fw);
            }
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void restoreObjects(List<Statefull> objects) {
        try {
            BufferedReader br = new BufferedReader(new FileReader("glass.state"));
            for (Statefull s : objects) {
                if (s!=null) s.restoreState(br);
            }
            br.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static class StructureFactorStuff implements Statefull {
        public final MeterStructureFactor meter;
        public final AccumulatorAverageFixed acc;
        public final DataSinkBlockAveragerSFac averager;
        public final StructureFactorComponentWriter writer;
        public StructureFactorStuff(MeterStructureFactor meter, AccumulatorAverageFixed acc, DataSinkBlockAveragerSFac averager, StructureFactorComponentWriter writer) {
            this.meter = meter;
            this.acc = acc;
            this.writer = writer;
            this.averager = averager;
        }

        @Override
        public void saveState(Writer fw) throws IOException {
            if (acc != null) acc.saveState(fw);
            if (writer != null) writer.saveState(fw);
            if (averager != null) averager.saveState(fw);
        }

        @Override
        public void restoreState(BufferedReader br) throws IOException {
            if (acc != null) acc.restoreState(br);
            if (writer != null) writer.restoreState(br);
            if (averager != null) averager.restoreState(br);
        }
    }
    public static class StructureFactorStuff2 implements Statefull {
        public final MeterStructureFactor meter;
        public final DataSinkBlockAveragerSFac averager;
        public final StructureFactorComponentCorrelation cor;
        public StructureFactorStuff2(MeterStructureFactor meter, DataSinkBlockAveragerSFac averager, StructureFactorComponentCorrelation cor) {
            this.meter = meter;
            this.averager = averager;
            this.cor = cor;
        }

        @Override
        public void saveState(Writer fw) throws IOException {
            if (averager != null) averager.saveState(fw);
            if (cor != null) cor.saveState(fw);
        }

        @Override
        public void restoreState(BufferedReader br) throws IOException {
            if (averager != null) averager.restoreState(br);
            if (cor != null) cor.restoreState(br);
        }
    }

    public static StructureFactorStuff2 setupStructureFactor(Box box, Vector[] wv, MeterStructureFactor.AtomSignalSource atomSignal,
                                                            int interval, ConfigurationStorage configStorage) {

        MeterStructureFactor meterSFac = new MeterStructureFactor(box, 1, atomSignal);
        meterSFac.setWaveVec(wv);
        meterSFac.setNormalizeByN(true);
        StructureFactorComponentCorrelation sfcDensityCor = new StructureFactorComponentCorrelation(2, configStorage);
        sfcDensityCor.setMinInterval(interval);
        DataSinkBlockAveragerSFac averager = new DataSinkBlockAveragerSFac(configStorage, interval, meterSFac);
        averager.addSink(sfcDensityCor);
        DataPump pumpSFac = new DataPump(meterSFac, averager);
        ConfigurationStoragePumper csp = new ConfigurationStoragePumper(pumpSFac, configStorage);
        configStorage.addListener(csp);
        csp.setPrevStep(interval);
        StructureFactorComponentCorrelation sfacCor = new StructureFactorComponentCorrelation(2, configStorage);
        averager.addSink(sfacCor);
        return new StructureFactorStuff2(meterSFac, averager, sfacCor);
    }

    public static StructureFactorStuff setupStructureFactor(Box box, double cut, MeterStructureFactor.AtomSignalSource atomSignal,
                                            int interval, ConfigurationStorage configStorage) {

        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);
        MeterStructureFactor meterSFac = new MeterStructureFactor(box, cut, atomSignal);
        meterSFac.setNormalizeByN(true);
        DataFork forkSFac = new DataFork();
        DataPump pumpSFac = new DataPump(meterSFac, forkSFac);
        ConfigurationStoragePumper csp = new ConfigurationStoragePumper(pumpSFac, configStorage);
        configStorage.addListener(csp);
        csp.setPrevStep(interval);
        forkSFac.addDataSink(accSFac);
        StructureFactorComponentWriter sfacWriter = null;
        DataSinkBlockAveragerSFac averager = new DataSinkBlockAveragerSFac(configStorage, interval, meterSFac);
        forkSFac.addDataSink(averager);
        sfacWriter = new StructureFactorComponentWriter(meterSFac, configStorage, interval, 500);
        averager.addSink(sfacWriter);
        return new StructureFactorStuff(meterSFac, accSFac, averager, sfacWriter);
    }

    public static class StructureFactorMobilityStuff {
        public final MeterStructureFactor meter;
        public final AccumulatorAverageFixed acc;
        public final StructureFactorComponentWriter writer;
        public StructureFactorMobilityStuff(MeterStructureFactor meter, AccumulatorAverageFixed acc, StructureFactorComponentWriter writer) {
            this.meter = meter;
            this.acc = acc;
            this.writer = writer;
        }
    }
    public static class StructureFactorMobilityStuff2 {
        public final MeterStructureFactor meter;
        public final DataFork fork;

        public StructureFactorMobilityStuff2(MeterStructureFactor meter, DataFork fork) {
            this.meter = meter;
            this.fork = fork;
        }
    }

    public static StructureFactorMobilityStuff setupMobilityStructureFactor(ConfigurationStorage configStorage, StructureFactorStuff sfacDensity, int N, AtomType otherType,
                                                                            int interval, double cut, Box box, StructureFactorComponentWriter writer) {
        AtomSignalMobility signalMobility = new AtomSignalMobility(configStorage, sfacDensity.meter, N);
        sfacDensity.averager.addSink(signalMobility);
        signalMobility.setAtomTypeFactor(otherType, 0);
        signalMobility.setPrevConfig(interval + 1);
        return setupMoMoStructureFactor(signalMobility, configStorage, interval, cut, box, writer);
    }

    public static StructureFactorMobilityStuff setupMotionStructureFactor(ConfigurationStorage configStorage, int xyz,
                                                                            int interval, double cut, Box box, StructureFactorComponentWriter writer) {
        AtomSignalMotion signalMobility = new AtomSignalMotion(configStorage, xyz);
        signalMobility.setPrevConfig(interval + 1);
        return setupMoMoStructureFactor(signalMobility, configStorage, interval, cut, box, writer);
    }

    public static StructureFactorMobilityStuff setupMoMoStructureFactor(MeterStructureFactor.AtomSignalSource signal, ConfigurationStorage configStorage,
                                                                        int interval, double cut, Box box, StructureFactorComponentWriter writer) {

        AtomPositionMobility positionMobility = new AtomPositionMobility(configStorage);
        positionMobility.setPrevConfig(interval + 1);
        MeterStructureFactor meter = new MeterStructureFactor(box, cut, signal, positionMobility);
        meter.setNormalizeByN(true);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1);
        DataFork fork = new DataFork();
        DataPump pumpSFacMobility = new DataPump(meter, fork);
        fork.addDataSink(acc);
        // ensures pump fires when config with delta t is available
        ConfigurationStoragePumper cspMobility = new ConfigurationStoragePumper(pumpSFacMobility, configStorage);
        cspMobility.setPrevStep(interval);
        configStorage.addListener(cspMobility);

        if (writer == null) writer = new StructureFactorComponentWriter(meter, configStorage, interval, 500);
        fork.addDataSink(new StructorFactorComponentExtractor(meter, interval, writer));

        return new StructureFactorMobilityStuff(meter, acc, writer);
    }


    public static StructureFactorMobilityStuff2 setupMobilityStructureFactor(ConfigurationStorage configStorage, StructureFactorStuff sfacDensity, int N, AtomType otherType,
                                                                            int interval, Vector[] wv, Box box, StructureFactorComponentCorrelation sfcCor) {
        AtomSignalMobility signalMobility = new AtomSignalMobility(configStorage, sfacDensity.meter, N);
        sfacDensity.averager.addSink(signalMobility);
        signalMobility.setAtomTypeFactor(otherType, 0);
        signalMobility.setPrevConfig(interval + 1);
        return setupMoMoStructureFactor(signalMobility, configStorage, interval, wv, box, sfcCor);
    }

    public static StructureFactorMobilityStuff2 setupMotionStructureFactor(ConfigurationStorage configStorage, int xyz,
                                                                          int interval, Vector[] wv, Box box, StructureFactorComponentCorrelation sfcCor) {
        AtomSignalMotion signalMobility = new AtomSignalMotion(configStorage, xyz);
        signalMobility.setPrevConfig(interval + 1);
        return setupMoMoStructureFactor(signalMobility, configStorage, interval, wv, box, sfcCor);
    }

    public static StructureFactorMobilityStuff2 setupMoMoStructureFactor(MeterStructureFactor.AtomSignalSource signal, ConfigurationStorage configStorage,
                                                                        int interval, Vector[] wv, Box box, StructureFactorComponentCorrelation sfcCor) {

        AtomPositionMobility positionMobility = new AtomPositionMobility(configStorage);
        positionMobility.setPrevConfig(interval + 1);
        MeterStructureFactor meter = new MeterStructureFactor(box, 1, signal, positionMobility);
        meter.setWaveVec(wv);
        meter.setNormalizeByN(true);
        DataFork fork = new DataFork();
        DataPump pumpSFacMobility = new DataPump(meter, fork);
        // ensures pump fires when config with delta t is available
        ConfigurationStoragePumper cspMobility = new ConfigurationStoragePumper(pumpSFacMobility, configStorage);
        cspMobility.setPrevStep(interval);
        configStorage.addListener(cspMobility);

        fork.addDataSink(sfcCor.makeSink(interval, meter));

        return new StructureFactorMobilityStuff2(meter, fork);
    }

    public static void main(String[] args) {
        double startTime = System.nanoTime()/1e9;
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.doLinear = !true;
            params.potential = SimGlass.PotentialChoice.HS;
            params.nA = 100;
            params.nB = 100;
            params.density =  1.5; // 2D params.density = 0.509733; //3D  = 0.69099;
            params.D = 3;
            params.temperature = 1.0;
            params.numStepsEqIsothermal = 100000;
            params.numStepsIsothermal = 0;
            params.numSteps = 100000;
            params.minDrFilter = 0.4;
            params.qx = new double[]{7.0};
            params.rcLJ = 2.5;
            params.randomSeeds = new int[]{-1329042323, 760258263, 1926332026, -524926192};
            params.tStep = 0.01;
        }

        double rho = params.density;
        if (params.eta>0 && params.potential == SimGlass.PotentialChoice.HS) {
            rho = SimGlass.densityForEta(params.eta, params.nA/(double)(params.nA+params.nB), params.sigmaB);
            params.density = rho;
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep, params.randomSeeds, params.rcLJ, params.sigmaB);
        if (params.randomSeeds == null) System.out.println("random seeds: " + Arrays.toString(sim.getRandomSeeds()));
        else System.out.println("set random seeds: " + Arrays.toString(params.randomSeeds));
        System.out.println(params.D + "D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        int numAtoms = params.nA + params.nB;
        System.out.println("T = " + params.temperature);
        System.out.println(params.numSteps + " production MD steps after " + params.numStepsEqIsothermal + " NVT equilibration steps, "+params.numStepsIsothermal+" NVT E check using dt = " + sim.integrator.getTimeStep());
        double volume = sim.box.getBoundary().volume();
        if (params.potential == SimGlass.PotentialChoice.HS) {
            double phi;
            double sb = params.sigmaB == 0 ? 1.0/1.4 : params.sigmaB;
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

        StepsState stepsState = new StepsState();
        StepsState savedSteps = new StepsState();
        int savedStage = 0; // last stage from previous run.  0 means no previous run
        // we need to know how far we got before reading glass.state so that we know what objects to try to read
        boolean haveSavedState = new File("glass.steps").exists();
        if (haveSavedState) {
            try {
                BufferedReader br = new BufferedReader(new FileReader("glass.steps"));
                savedSteps.restoreState(br);
                br.close();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            if (savedSteps.numStepsEqIsothermal > 0) savedStage = 1; // eq ran
            if (savedSteps.numStepsIsothermal > 0) savedStage = 2;   // nvt ran
            if (savedSteps.numSteps > 0) savedStage = 3;             // production ran
            stepsState.numStepsEqIsothermal = savedSteps.numStepsEqIsothermal;
            stepsState.numStepsIsothermal = savedSteps.numStepsIsothermal;
            stepsState.numSteps = savedSteps.numSteps;
            params.numStepsEqIsothermal -= savedSteps.numStepsEqIsothermal;
            params.numStepsIsothermal -= savedSteps.numStepsIsothermal;
            params.numSteps -= savedSteps.numSteps;
        }
        List<Statefull> objects = new ArrayList<>();
        objects.add(sim.box);
        objects.add(sim.integrator);
        boolean skipReset = false;
        if (params.numStepsEqIsothermal > 0 && haveSavedState && savedStage == 1) {
            // we have a restore file and we need to do isothermal equilibration
            // we should restore and then continue equilibrating isothermally
            restoreObjects(objects);
            System.out.println("Continuing equilibration after " + sim.integrator.getStepCount() + " steps");
            // find neighbors with new config
            // reset collision times (for hard) or compute forces (for soft)
            sim.integrator.postRestore();
            skipReset = true; // allow integrator to retain its step count
        }

        //Equilibration
        long time0eq = System.nanoTime();
        double temperature0 = Math.max(params.temperatureMelt, params.temperature);
        if (temperature0 > params.temperature) System.out.println("Equilibrating at T=" + temperature0);
        sim.integrator.setIsothermal(true);
        if (params.doSwap) sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(temperature0);
        // Equilibrate isothermally
        double maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsEqIsothermal, false, maxTime, skipReset));
        if (System.nanoTime()/1e9 - startTime >= params.maxWalltime) {
            // we've reached our time limit.  disable the rest of the sim.
            params.numSteps = params.numStepsIsothermal = 0;
            // and set our eq steps to what we were able to actually run
            System.out.println(sim.integrator.getStepCount()+" steps of isothermal equilibration finished before time limit");
        }

        long time1eq = System.nanoTime();
        // if we actually ran eq, then update the with # of steps we've finished
        if (params.numStepsEqIsothermal>0) stepsState.numStepsEqIsothermal = sim.integrator.getStepCount();

        if (params.numStepsIsothermal + params.numSteps == 0) {
            // we're done!
            saveObjects(stepsState, objects);

            System.out.printf("\nisothermal equilibration time: %3.2f s\n", (time1eq-time0eq)/1e9);
            return;
        }

        skipReset = false;
        if (params.numStepsEqIsothermal == 0 && haveSavedState && savedStage == 1) {
            // we have a restore file, we did not do any isothermal equilibration, but we need to do isothermal
            // collection of total energy
            // we should restore and then collect total energy
            // we did not start collecting E, so we need to restore before we add accE to our objects
            restoreObjects(objects);
            System.out.println("Running NVT after " + sim.integrator.getStepCount() + " equilibration steps");
            // find neighbors with new config
            // reset collision times (for hard) or compute forces (for soft)
            sim.integrator.postRestore();
        }

        AccumulatorAverageFixed accE = new AccumulatorAverageFixed(1);
        MeterEnergyFromIntegrator meterE = new MeterEnergyFromIntegrator(sim.integrator);
        DataPumpListener pumpE = new DataPumpListener(meterE, accE, 10);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            sim.integrator.getEventManager().addListener(pumpE);
            // add now to read from glass.state
            objects.add(accE);
        }
        if (params.numStepsEqIsothermal == 0 && params.numStepsIsothermal > 0 && haveSavedState && savedStage == 2) {
            // we have a restore file, we previously did isothermal collection of total E and need to continue that
            // we should restore and then collect total energy
            restoreObjects(objects);
            System.out.println("Continuing NVT stage after " + sim.integrator.getStepCount() + " steps");
            // find neighbors with new config
            // reset collision times (for hard) or compute forces (for soft)
            sim.integrator.postRestore();
            skipReset = true;
        }

        if (temperature0 == params.temperature) {
            // Now collect average total energy (stage 2).
            maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsIsothermal, false, maxTime, skipReset));
            if (System.nanoTime()/1e9 - startTime >= params.maxWalltime) {
                // we've reached our time limit.  disable the rest of the sim.
                params.numSteps = 0;
                // and set our nvt steps to what we were able to actually run
                System.out.println(sim.integrator.getStepCount()+" steps of isothermal finished before time limit");
            }
            if (params.numStepsIsothermal>0) stepsState.numStepsIsothermal = sim.integrator.getStepCount();
        }
        else {
            // Don't be trying to save/restore into here.
            // Just do a run at a high T, and then another at the T you want.

            // first stage was at higher temperature to ensure crystal melted
            // now we need a second stage at the temperature we actually want
            System.out.println("Equilibrating at T=" + params.temperature);
            sim.integrator.setTemperature(params.temperature);
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsEqIsothermal));
            stepsState.numStepsEqIsothermal += params.numStepsEqIsothermal;
            accE.reset();

            // Assume isothermal equilibration.  Now collect average total energy.
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numStepsIsothermal));
            stepsState.numStepsEqIsothermal += params.numStepsIsothermal;
        }

        long time2eq = System.nanoTime();
        if (params.numStepsEqIsothermal + params.numStepsIsothermal > 0 && params.doSwap) {
            // only for this sim, swapMove not saved/restored (has no state, only tracker has state)
            System.out.println("swap acceptance: "+sim.swapMove.getTracker().acceptanceProbability());
        }
        System.out.printf("\nequilibration time: %3.2f s\n", (time2eq-time0eq)/1e9);
        if (params.numSteps == 0) {
            saveObjects(stepsState, objects);

            return;
        }

        if (params.numStepsEqIsothermal == 0 && params.numStepsIsothermal == 0 && haveSavedState && savedStage < 3) {
            // we have a restore file, we need to do production, but have never done production before.
            // we need to restore here before setting up our data collection
            restoreObjects(objects);
            System.out.println("Running production after continuing from previous eq sim");
            // find neighbors with new config
            // reset collision times (for hard) or compute forces (for soft)
            sim.integrator.postRestore();
        }

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

        //Linear Viscosity
        AccumulatorLinearPTensor pTensorLinearAccumVisc = null;
        if (params.doLinear) {
            pTensorLinearAccumVisc = new AccumulatorLinearPTensor(sim.box, dn * sim.integrator.getTimeStep(), 1000);
            pTensorLinearAccumVisc.setEnabled(true);
            pTensorFork.addDataSink(pTensorLinearAccumVisc);
        }

        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork, dn);


        AccumulatorAverageFixed gTensorAccumulator = null;
        if (params.potential != SimGlass.PotentialChoice.HS) {
            gTensorAccumulator = new AccumulatorAverageFixed(blocksize);
            DataProcessor dpSquared = new DataProcessor() {
                final DataDouble data = new DataDouble();

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
            dpSquared.setDataSink(gTensorAccumulator);
        }

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
            tAccumulator.addDataSink(pTensorLinearAccumVisc.makeTemperatureSink(), new AccumulatorAverage.StatType[]{tAccumulator.AVERAGE});
            sim.integrator.getEventManager().addListener(tPump);
        }


        //shear stress AC
        AccumulatorAutocorrelationShearStress dpxyAutocor = new AccumulatorAutocorrelationShearStress(256, sim.integrator.getTimeStep());
        boolean doPxyAutocor = params.doPxyAutocor && sim.potentialChoice != SimGlass.PotentialChoice.HS;
        if (doPxyAutocor) {
            pTensorFork.addDataSink(dpxyAutocor);
        }
        dpxyAutocor.setPushInterval(16384);

        //linear MSD
        ConfigurationStorage configStorageLinearMSD = null;
        DataSourceLinearMSD meterLinearMSD = null;
        if (params.doLinear){
            configStorageLinearMSD = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.LINEAR, 1000, 1);
            configStorageLinearMSD.setEnabled(true);
            sim.integrator.getEventManager().addListener(configStorageLinearMSD);
            meterLinearMSD = new DataSourceLinearMSD(configStorageLinearMSD);
            configStorageLinearMSD.addListener(meterLinearMSD);
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

        DataSourceCorBlock dsCorMSD = new DataSourceCorBlock(sim.integrator);
        dsCorMSD.setMinInterval(3);
        meterMSD.addMSDSink(dsCorMSD);
        dsCorMSD.setEnabled(true);
        DataSourceBlockAvgCor dsCorP = new DataSourceBlockAvgCor(sim.integrator);
        pFork.addDataSink(dsCorP);
        dsCorP.setEnabled(true);
        DataSourceMSDcorP dsMSDcorP = new DataSourceMSDcorP(sim.integrator);
        pFork.addDataSink(dsMSDcorP);
        meterMSD.addMSDSink(dsMSDcorP);
        dsMSDcorP.setEnabled(true);

        DataSourceCorBlock dsCorVisc = new DataSourceCorBlock(sim.integrator);
        pTensorAccumVisc.addViscositySink(dsCorVisc);
        dsCorVisc.setEnabled(true);

        //Linear-VAC
        DataSourceLinearVAC meterLinearVAC = null;
        if (params.doLinear){
            configStorageLinearMSD.setDoVelocity(true);
            meterLinearVAC = new DataSourceLinearVAC(configStorageLinearMSD);
            configStorageLinearMSD.addListener(meterLinearVAC);
        }

        //VAC
        configStorageMSD.setDoVelocity(true);
        DataSourceVAC meterVAC = new DataSourceVAC(configStorageMSD);
        configStorageMSD.addListener(meterVAC);

        //Fs
        DataSourceFs[][] meterFs = new DataSourceFs[params.qx.length][2];
        DataSourceF[] meterF = new DataSourceF[params.qx.length];
        for (int i=0; i<meterFs.length; i++) {
            Vector q = sim.getSpace().makeVector();
            q.setX(0, params.qx[i]);
            for (int j=0; j<2; j ++) {
                meterFs[i][j] = new DataSourceFs(configStorageMSD);
                meterFs[i][j].setQ(q);
                meterFs[i][j].setAtomType(j==0 ? sim.speciesA.getLeafType() : sim.speciesB.getLeafType());
                configStorageMSD.addListener(meterFs[i][j]);
            }

            meterF[i] = new DataSourceF(configStorageMSD, params. qx[i]);
            configStorageMSD.addListener(meterF[i]);
        }

        //Percolation
        int percMinInterval = 11;
        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorageMSD);
        atomFilterDeviation.setMinDistance(params.minDrFilter);
        atomFilterDeviation.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(sim.getSpeciesManager(), configStorageMSD, atomFilterDeviation, percMinInterval);
        configStorageMSD.addListener(meterPerc);

        AtomTestDeviation atomFilterDeviation3 = new AtomTestDeviation(sim.box, configStorageMSD3);
        atomFilterDeviation3.setMinDistance(params.minDrFilter);
        atomFilterDeviation3.setDoMobileOnly(false);
        DataSourcePercolation meterPerc3 = new DataSourcePercolation(sim.getSpeciesManager(), configStorageMSD3, atomFilterDeviation3, percMinInterval - 1, meterPerc.getHistogram());
        configStorageMSD3.addListener(meterPerc3);

        DataSourcePercolation0 meterPerc0 = new DataSourcePercolation0(sim.box, sim.getRandom());
        meterPerc0.setImmFracs(new double[]{0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27});
        AccumulatorAverageFixed accPerc0 = new AccumulatorAverageFixed(10);
        DataPumpListener pumpPerc0 = new DataPumpListener(meterPerc0, accPerc0, 10000);
        sim.integrator.getEventManager().addListener(pumpPerc0);

        DataSourceQ4 meterQ4 = new DataSourceQ4(configStorageMSD, percMinInterval + 1);
        meterQ4.setMaxDr(params.minDrFilter);
        configStorageMSD.addListener(meterQ4);

        //Strings
        DataSourceStrings meterL = new DataSourceStrings(configStorageMSD, 6);
        configStorageMSD.addListener(meterL);

        //Alpha2
        DataSourceAlpha2 meterAlpha2A = new DataSourceAlpha2(configStorageMSD, sim.speciesA);
        configStorageMSD.addListener(meterAlpha2A);
        DataSourceAlpha2 meterAlpha2B = new DataSourceAlpha2(configStorageMSD, sim.speciesB);
        configStorageMSD.addListener(meterAlpha2B);

        //structure factors for all low-wavelength WV.  we just collect simple averages here
        double cut3 = 3;
        if (numAtoms > 800) cut3 /= Math.pow(numAtoms / 500.0, 1.0 / sim.getSpace().D());
        MeterStructureFactor.AtomSignalSource atomSignalSimple = new MeterStructureFactor.AtomSignalSourceByType();

        StructureFactorStuff sfacDensity = setupStructureFactor(sim.box, cut3, atomSignalSimple, params.sfacMinInterval, configStorageMSD);

        StructureFactorStuff sfacA = sfacDensity, sfacB = null;
        StructureFactorStuff sfacPack = null, sfacAB = null;
        StructureFactorStuff sfacStress = null;

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
        AtomSignalStress signalStressNormal = new AtomSignalStress(stressSource, normalComps);

        if (params.nB>0) {
            MeterStructureFactor.AtomSignalSourceByType atomSignalA = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalA.setAtomTypeFactor(sim.speciesB.getLeafType(), 0);
            sfacA = setupStructureFactor(sim.box, cut3, atomSignalA, params.sfacMinInterval, configStorageMSD);

            MeterStructureFactor.AtomSignalSourceByType atomSignalB = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalB.setAtomTypeFactor(sim.speciesA.getLeafType(), 0);
            sfacB = setupStructureFactor(sim.box, cut3, atomSignalB, params.sfacMinInterval, configStorageMSD);

            double vA = sim.getSpace().sphereVolume(0.5);
            double vB = sim.getSpace().sphereVolume(0.5*sim.sigmaB);

            MeterStructureFactor.AtomSignalSourceByType atomSignalPack = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalPack.setAtomTypeFactor(sim.speciesA.getAtomType(0), +vA);
            atomSignalPack.setAtomTypeFactor(sim.speciesB.getAtomType(0), +vB);
            sfacPack = setupStructureFactor(sim.box, cut3, atomSignalPack, params.sfacMinInterval, configStorageMSD);

            MeterStructureFactor.AtomSignalSourceByType atomSignalAB = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalAB.setAtomTypeFactor(sim.speciesA.getAtomType(0), +vA);
            atomSignalAB.setAtomTypeFactor(sim.speciesB.getAtomType(0), -vB);
            sfacAB = setupStructureFactor(sim.box, cut3, atomSignalAB, params.sfacMinInterval, configStorageMSD);

            sfacStress = setupStructureFactor(sim.box, cut3, signalStressNormal, params.sfacMinInterval, configStorageMSD);
        }

        AtomSignalKineticEnergy atomSignalKE = new AtomSignalKineticEnergy();
        StructureFactorStuff sfacKE = setupStructureFactor(sim.box, cut3, atomSignalKE, params.sfacMinInterval, configStorageMSD);

        StructureFactorMobilityStuff[] sfacMobilityA = new StructureFactorMobilityStuff[30];
        StructureFactorMobilityStuff[] sfacMobilityB = new StructureFactorMobilityStuff[30];
        StructureFactorMobilityStuff[] sfacMotionX = new StructureFactorMobilityStuff[30];
        StructureFactorComponentWriter sfacMobilityWriterA = null, sfacMobilityWriterB = null, sfacMotionXWriter = null;
        for (int i = params.sfacMinInterval; i < 30; i++) {
            sfacMobilityA[i] = setupMobilityStructureFactor(configStorageMSD, sfacA, params.nA, sim.speciesB.getLeafType(), i, cut3, sim.box, sfacMobilityWriterA);
            sfacMobilityWriterA = sfacMobilityA[i].writer;

            if (params.nB>0) {
                sfacMobilityB[i] = setupMobilityStructureFactor(configStorageMSD, sfacB, params.nB, sim.speciesA.getLeafType(), i, cut3, sim.box, sfacMobilityWriterB);
                sfacMobilityWriterB = sfacMobilityB[i].writer;
            }

            sfacMotionX[i] = setupMotionStructureFactor(configStorageMSD, 0, i, cut3, sim.box, sfacMotionXWriter);
            sfacMotionXWriter = sfacMotionX[i].writer;
        }

        // now rinse and repeat with only 2 wave vectors.  we will use these to compute correlations
        // for motion, this is a compression mode
        Vector[] wvx = new Vector[2];
        double wv1 = 2*Math.PI/sim.box.getBoundary().getBoxSize().getX(0);
        wvx[0] = sim.getSpace().makeVector();
        wvx[1] = sim.getSpace().makeVector();
        wvx[0].setX(0, wv1);
        wvx[1].setX(0, 2*wv1);

        // only need this for motion, this is a shear mode
        Vector[] wvy = new Vector[2];
        wvy[0] = sim.getSpace().makeVector();
        wvy[1] = sim.getSpace().makeVector();
        wvy[0].setX(1, wv1);
        wvy[1].setX(1, 2*wv1);

        StructureFactorStuff2 sfacDensityX = setupStructureFactor(sim.box, wvx, atomSignalSimple, 1, configStorageMSD);

        StructureFactorStuff2 sfacAX = sfacDensityX, sfacBX = null;
        StructureFactorStuff2 sfacPackX = null, sfacABX = null;
        StructureFactorStuff2 sfacStressX = null;

        if (params.nB>0) {
            MeterStructureFactor.AtomSignalSourceByType atomSignalA = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalA.setAtomTypeFactor(sim.speciesB.getLeafType(), 0);
            sfacAX = setupStructureFactor(sim.box, wvx, atomSignalA, params.sfacMinInterval, configStorageMSD);

            MeterStructureFactor.AtomSignalSourceByType atomSignalB = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalB.setAtomTypeFactor(sim.speciesA.getLeafType(), 0);
            sfacBX = setupStructureFactor(sim.box, wvx, atomSignalB, params.sfacMinInterval, configStorageMSD);

            double vA = sim.getSpace().sphereVolume(0.5);
            double vB = sim.getSpace().sphereVolume(0.5*sim.sigmaB);

            MeterStructureFactor.AtomSignalSourceByType atomSignalPack = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalPack.setAtomTypeFactor(sim.speciesA.getAtomType(0), +vA);
            atomSignalPack.setAtomTypeFactor(sim.speciesB.getAtomType(0), +vB);
            sfacPackX = setupStructureFactor(sim.box, wvx, atomSignalPack, params.sfacMinInterval, configStorageMSD);

            MeterStructureFactor.AtomSignalSourceByType atomSignalAB = new MeterStructureFactor.AtomSignalSourceByType();
            atomSignalAB.setAtomTypeFactor(sim.speciesA.getAtomType(0), +vA);
            atomSignalAB.setAtomTypeFactor(sim.speciesB.getAtomType(0), -vB);
            sfacABX = setupStructureFactor(sim.box, wvx, atomSignalAB, params.sfacMinInterval, configStorageMSD);

            sfacStressX = setupStructureFactor(sim.box, wvx, signalStressNormal, params.sfacMinInterval, configStorageMSD);
        }

        StructureFactorStuff2 sfacKEX = setupStructureFactor(sim.box, wvx, atomSignalKE, params.sfacMinInterval, configStorageMSD);

        StructureFactorMobilityStuff2[] sfacMobilityAX = new StructureFactorMobilityStuff2[30];
        StructureFactorMobilityStuff2[] sfacMobilityBX = new StructureFactorMobilityStuff2[30];
        StructureFactorMobilityStuff2[] sfacMotionXX = new StructureFactorMobilityStuff2[30];
        StructureFactorMobilityStuff2[] sfacMotionXY = new StructureFactorMobilityStuff2[30];
        StructureFactorComponentCorrelation sfcMobilityACor = new StructureFactorComponentCorrelation(2, configStorageMSD);
        sfcMobilityACor.setMinInterval(params.sfacMinInterval);
        StructureFactorComponentCorrelation sfcMobilityBCor = new StructureFactorComponentCorrelation(2, configStorageMSD);
        sfcMobilityBCor.setMinInterval(params.sfacMinInterval);
        StructureFactorComponentCorrelation sfcMotionXXCor = new StructureFactorComponentCorrelation(2, configStorageMSD);
        sfcMotionXXCor.setMinInterval(params.sfacMinInterval);
        StructureFactorComponentCorrelation sfcMotionXYCor = new StructureFactorComponentCorrelation(2, configStorageMSD);
        sfcMotionXYCor.setMinInterval(params.sfacMinInterval);
        for (int i = params.sfacMinInterval; i < 30; i++) {
            sfacMobilityAX[i] = setupMobilityStructureFactor(configStorageMSD, sfacA, params.nA, sim.speciesB.getLeafType(), i, wvx, sim.box, sfcMobilityACor);

            if (params.nB>0) {
                sfacMobilityBX[i] = setupMobilityStructureFactor(configStorageMSD, sfacB, params.nB, sim.speciesA.getLeafType(), i, wvx, sim.box, sfcMobilityBCor);
            }

            sfacMotionXX[i] = setupMotionStructureFactor(configStorageMSD, 0, i, wvx, sim.box, sfcMotionXXCor);
            sfacMotionXY[i] = setupMotionStructureFactor(configStorageMSD, 0, i, wvx, sim.box, sfcMotionXYCor);
        }


        DataSourceCorrelation dsCorSFacPackMobilityA = new DataSourceCorrelation(configStorageMSD, 2);
        DataSourceCorrelation dsCorSFacPackMobilityB = null, dsCorSFacABMobilityA = null, dsCorSFacABMobilityB = null;
        DataSourceCorrelation dsCorSFacStressMobilityA = new DataSourceCorrelation(configStorageMSD, 2);
        DataSourceCorrelation dsCorSFacStressMobilityB = null;
        DataSourceCorrelation dsCorSFacPackPackAB = null;
        DataSourceCorrelation dsCorSFacKEMobilityB = null;
        if (params.nB>0) {
            sfacPackX.averager.addSink(dsCorSFacPackMobilityA.makeReceiver(0));
            dsCorSFacPackMobilityB = new DataSourceCorrelation(configStorageMSD, 2);
            sfacPackX.averager.addSink(dsCorSFacPackMobilityB.makeReceiver(0));
            dsCorSFacABMobilityA = new DataSourceCorrelation(configStorageMSD, 2);
            sfacABX.averager.addSink(dsCorSFacABMobilityA.makeReceiver(0));
            dsCorSFacABMobilityB = new DataSourceCorrelation(configStorageMSD, 2);
            sfacABX.averager.addSink(dsCorSFacABMobilityB.makeReceiver(0));
            dsCorSFacPackPackAB = new DataSourceCorrelation(configStorageMSD, 2);
            dsCorSFacStressMobilityB = new DataSourceCorrelation(configStorageMSD, 2);
            sfacStressX.averager.addSink(dsCorSFacStressMobilityA.makeReceiver(0));
            sfacStressX.averager.addSink(dsCorSFacStressMobilityB.makeReceiver(0));
            sfacPackX.averager.addSink(dsCorSFacPackPackAB.makeReceiver(0));
            sfacABX.averager.addSink(dsCorSFacPackPackAB.makeReceiver(1));
            dsCorSFacKEMobilityB = new DataSourceCorrelation(configStorageMSD, 2);
            sfacKEX.averager.addSink(dsCorSFacKEMobilityB.makeReceiver(0));
        }
        else {
            // pack and density are equivalent
            sfacDensityX.averager.addSink(dsCorSFacPackMobilityA.makeReceiver(0));
            sfacStressX.averager.addSink(dsCorSFacPackMobilityA.makeReceiver(0));
        }


        DataSourceCorrelation dsCorSFacKEMobilityA = new DataSourceCorrelation(configStorageMSD, 2);
        sfacKEX.averager.addSink(dsCorSFacKEMobilityA.makeReceiver(0));

        DataSourceCorrelation dsCorSFacMobilityAB = new DataSourceCorrelation(configStorageMSD, 2);

        for (int i = params.sfacMinInterval; i < 30; i++) {
            sfacMobilityAX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityAX[i].meter, i, dsCorSFacPackMobilityA));
            sfacMobilityAX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityAX[i].meter, i, dsCorSFacABMobilityA));
            sfacMobilityAX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityAX[i].meter, i, dsCorSFacKEMobilityA));
            sfacMobilityAX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityAX[i].meter, i, dsCorSFacStressMobilityA));

            if (params.nB>0) {
                sfacMobilityBX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityBX[i].meter, i, dsCorSFacPackMobilityB));
                sfacMobilityBX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityBX[i].meter, i, dsCorSFacABMobilityB));
                sfacMobilityBX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityBX[i].meter, i, dsCorSFacKEMobilityB));
                sfacMobilityBX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityBX[i].meter, i, dsCorSFacStressMobilityB));

                sfacMobilityAX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityAX[i].meter, i, dsCorSFacMobilityAB, 0));
                sfacMobilityBX[i].fork.addDataSink(new StructorFactorComponentExtractor(sfacMobilityBX[i].meter, i, dsCorSFacMobilityAB));
            }
        }

        double xGsMax = 3;
        int gsMinConfig = 6;
        MeterGs meterGs = new MeterGs(configStorageMSD);
        meterGs.setMinConfigIndex(gsMinConfig);
        configStorageMSD.addListener(meterGs);
        meterGs.getXDataSource().setXMax(xGsMax);
        MeterGs meterGsA = new MeterGs(configStorageMSD);
        meterGsA.setAtomTypes(sim.speciesA.getLeafType());
        meterGsA.setMinConfigIndex(gsMinConfig);
        configStorageMSD.addListener(meterGsA);
        meterGsA.getXDataSource().setXMax(xGsMax);
        MeterGs meterGsB = new MeterGs(configStorageMSD);
        meterGsB.setAtomTypes(sim.speciesB.getLeafType());
        meterGsB.setMinConfigIndex(gsMinConfig);
        if (params.nB > 0) configStorageMSD.addListener(meterGsB);
        meterGsB.getXDataSource().setXMax(xGsMax);

        MeterCorrelationSelf meterCorrelationSelf = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.TOTAL);
        configStorageMSD.addListener(meterCorrelationSelf);
        MeterCorrelationSelf meterCorrelationSelfMagA = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagA.setAtomType(sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterCorrelationSelfMagA);
        MeterCorrelationSelf meterCorrelationSelfMagB = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagB.setAtomType(sim.speciesB.getLeafType());
        if (params.nB>0) configStorageMSD.addListener(meterCorrelationSelfMagB);

        CorrelationSelf2 correlationSelf2 = new CorrelationSelf2(configStorageMSD, CorrelationSelf2.CorrelationType.TOTAL, 0.001, 20);
        configStorageMSD.addListener(correlationSelf2);

        sim.integrator.getEventManager().addListener(configStorageMSD3);
        sim.integrator.getEventManager().addListener(configStorageMSD);

        objects.add(pTensorAccumVisc);
        if (gTensorAccumulator != null) {
            objects.add(gTensorAccumulator);
        }
        objects.add(pAccumulator);
        if (tAccumulator != null) {
            objects.add(tAccumulator);
            objects.add(accPE);
        }
        objects.add(configStorageMSD);
        objects.add(configStorageMSD3);
        objects.add(meterMSD);
        objects.add(meterMSDA);
        if (params.nB>0) objects.add(meterMSDB);
        objects.add(dsCorMSD);
        objects.add(dsCorVisc);
        objects.add(dsCorP);
        objects.add(dsMSDcorP);
        objects.add(meterVAC);
        if (params.doLinear){
            objects.add(configStorageLinearMSD);
            objects.add(meterLinearMSD);
            objects.add(meterLinearVAC);
            objects.add(pTensorLinearAccumVisc);
        }
        for (int i=0; i<meterFs.length; i++) {
            for (int j=0; j<2; j++) objects.add(meterFs[i][j]);
            objects.add(meterF[i]);
        }
        objects.add(meterPerc);
        objects.add(meterPerc3);
        objects.add(accPerc0);
        objects.add(meterQ4);
        objects.add(meterL);
        objects.add(meterAlpha2A);
        objects.add(meterAlpha2B);
        objects.add(sfacDensity.acc);
        objects.add(sfacA.acc);
        if (params.nB>0) objects.add(sfacB.acc);
        objects.add(sfacKE.acc);
        if (params.nB>0) objects.add(sfacAB.acc);
        for (int i = params.sfacMinInterval; i < 30; i++) {
            objects.add(sfacMobilityA[i].acc);
            if (params.nB>0) objects.add(sfacMobilityB[i].acc);
            objects.add(sfacMotionX[i].acc);
        }
        objects.add(sfacMobilityWriterA);
        if (params.nB>0) objects.add(sfacMobilityWriterB);
        objects.add(sfacMotionXWriter);

        objects.add(sfacDensity);
        objects.add(sfacA);
        objects.add(sfacB);
        objects.add(sfacPack);
        objects.add(sfacAB);
        objects.add(sfacKE);
        objects.add(sfacStress);
        objects.add(sfacDensityX);
        objects.add(sfacAX);
        objects.add(sfacBX);
        objects.add(sfacPackX);
        objects.add(sfacABX);
        objects.add(sfacKEX);
        objects.add(sfacStressX);

        objects.add(sfcMobilityACor);
        objects.add(sfcMobilityBCor);
        objects.add(sfcMotionXXCor);
        objects.add(sfcMotionXYCor);

        objects.add(dsCorSFacPackPackAB);
        objects.add(dsCorSFacPackMobilityA);
        objects.add(dsCorSFacPackMobilityB);
        objects.add(dsCorSFacABMobilityA);
        objects.add(dsCorSFacABMobilityB);
        objects.add(dsCorSFacKEMobilityA);
        objects.add(dsCorSFacKEMobilityB);
        objects.add(dsCorSFacStressMobilityA);
        objects.add(dsCorSFacStressMobilityB);
        objects.add(dsCorSFacMobilityAB);

        objects.add(meterGs);
        objects.add(meterGsA);
        if (params.nB>0) objects.add(meterGsB);
        objects.add(meterCorrelationSelf);
        objects.add(meterCorrelationSelfMagA);
        if (params.nB>0) objects.add(meterCorrelationSelfMagB);
        objects.add(correlationSelf2);

        skipReset = false;
        if (params.numStepsEqIsothermal == 0 && params.numStepsIsothermal == 0 && haveSavedState && savedStage == 3) {
            // we started a production run previously
            restoreObjects(objects);
            System.out.println("Continuing production after " + sim.integrator.getStepCount() + " steps");
            // find neighbors with new config
            // reset collision times (for hard) or compute forces (for soft)
            sim.integrator.postRestore();
            skipReset = true;
        }

        //Run
        long time0 = System.nanoTime();
        maxTime = params.maxWalltime - (System.nanoTime()/1e9 - startTime);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps, false, maxTime, skipReset));
        if (System.nanoTime()/1e9 - startTime >= params.maxWalltime) {
            // we've reached our time limit
            // set our production steps to what we were able to actually run
            params.numSteps = sim.integrator.getStepCount();
            System.out.println(params.numSteps+" steps of production finished before time limit");
        }
        if (params.numSteps>0) stepsState.numSteps = sim.integrator.getStepCount();
        saveObjects(stepsState, objects);

//        sfacDensity.writer.writeFile("sfacDensityTraj.dat");
        sfacAB.writer.writeFile("sfacABTraj.dat");
        sfacPack.writer.writeFile("sfacPackTraj.dat");
//        sfacKE.writer.writeFile("sfacKETraj.dat");
        sfacStress.writer.writeFile("sfacStressTraj.dat");
        sfacMobilityWriterA.writeFile("sfacMobilityATraj.dat");
        sfacMobilityWriterB.writeFile("sfacMobilityBTraj.dat");
        sfacMotionXWriter.writeFile("sfacMotionTraj.dat");

        //Pressure
        DataGroup dataP = (DataGroup) pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg = dataPAvg.getValue(0);
        double pErr = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);

        String filenameLinearVisc=null, filenameLinearMSD=null, filenameLinearVAC=null;
        String filenameVisc, filenameMSD, filenameMSDA, filenameMSDB, filenamePerc,
                filenamePerc0, filenameImmFrac, filenameImmFracA, filenameImmFracB, filenameImmFracPerc, filenameL, filenameAlpha2A, filenameAlpha2B, filenameVAC;
        String[][] filenameFs = new String[meterFs.length][2];
        String[] filenameF = new String[meterFs.length];

        String fileTag;
        String filenamePxyAC = "";
        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            System.out.println("Z: " + pAvg / params.density / params.temperature + "  " + pErr / params.density / params.temperature + "  cor: " + pCorr);
            fileTag = String.format("Rho%1.3f", rho);
            filenameVisc = String.format("viscRho%1.3f.dat", rho);
            if (params.doLinear){
                filenameLinearVisc = String.format("viscLinearRho%1.3f.dat", rho);
                filenameLinearMSD = String.format("msdLinearRho%1.3f.dat", rho);
                filenameLinearVAC = String.format("vacLinearRho%1.3f.dat", rho);
            }
            filenameMSD = String.format("msdRho%1.3f.dat", rho);
            filenameMSDA = String.format("msdARho%1.3f.dat", rho);
            filenameMSDB = String.format("msdBRho%1.3f.dat", rho);
            filenameVAC = String.format("vacRho%1.3f.dat", rho);

            for (int i=0; i<meterFs.length; i++) {
                filenameFs[i][0] = String.format("FsARho%1.3fQ%1.2f.dat", rho, params.qx[i]);
                filenameFs[i][1] = String.format("FsBRho%1.3fQ%1.2f.dat", rho, params.qx[i]);
                filenameF[i] = String.format("FRho%1.3fQ%1.2f.dat", rho, params.qx[i]);
            }
            filenamePerc = String.format("percRho%1.3f.dat", rho);
            filenamePerc0 = String.format("perc0Rho%1.3f.dat", rho);
            filenameL = String.format("lRho%1.3f.dat", rho);
            filenameImmFrac = String.format("immFracRho%1.3f.dat", rho);
            filenameImmFracA = String.format("immFracARho%1.3f.dat", rho);
            filenameImmFracB = String.format("immFracBRho%1.3f.dat", rho);
            filenameImmFracPerc = String.format("immFracPercRho%1.3f.dat", rho);
            filenameAlpha2A = String.format("alpha2ARho%1.3f.dat", rho);
            filenameAlpha2B = String.format("alpha2BRho%1.3f.dat", rho);
        } else {
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
            fileTag = String.format("Rho%1.3fT%1.3f", rho, params.temperature);
            filenameVisc = String.format("viscRho%1.3fT%1.3f.dat", rho, params.temperature);
            if (params.doLinear) {
                filenameLinearVisc = String.format("viscLinearRho%1.3fT%1.3f.dat", rho, params.temperature);
                filenameLinearMSD = String.format("msdLinearRho%1.3fT%1.3f.dat", rho, params.temperature);
                filenameLinearVAC = String.format("vacLinearRho%1.3fT%1.3f.dat", rho, params.temperature);
            }

            filenameMSD = String.format("msdRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameMSDA = String.format("msdARho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameMSDB = String.format("msdBRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameVAC = String.format("vacRho%1.3fT%1.3f.dat", rho, params.temperature);
            for (int i=0; i<meterFs.length; i++) {
                filenameFs[i][0] = String.format("FsARho%1.3fT%1.3fQ%1.2f.dat", rho, params.temperature, params.qx[i]);
                filenameFs[i][1] = String.format("FsBRho%1.3fT%1.3fQ%1.2f.dat", rho, params.temperature, params.qx[i]);
                filenameF[i] = String.format("FRho%1.3fQ%1.2f.dat", rho, params.qx[i]);
            }
            filenamePerc = String.format("percRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenamePerc0 = String.format("perc0Rho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameL = String.format("lRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameImmFrac = String.format("immFracRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameImmFracA = String.format("immFracARho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameImmFracB = String.format("immFracBRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameImmFracPerc = String.format("immFracPercRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameAlpha2A = String.format("alpha2ARho%1.3fT%1.3f.dat", rho, params.temperature);
            filenameAlpha2B = String.format("alpha2BRho%1.3fT%1.3f.dat", rho, params.temperature);
            filenamePxyAC = String.format("acPxyRho%1.3fT%1.3f.dat", rho, params.temperature);
            double V = sim.box.getBoundary().volume();
            System.out.println("G: " + V * avgG / tAvg + " " + V * errG / tAvg + " cor: " + corG + "\n");
        }
        System.out.println("P: " + pAvg + "  " + pErr + "  cor: " + pCorr);

        try {
            if (doPxyAutocor) {
                GlassProd.writeDataToFile(dpxyAutocor, filenamePxyAC);
            }
            GlassProd.writeDataToFile(dsCorMSD, "MSDcor.dat");
            GlassProd.writeDataToFile(dsCorP, "Pcor.dat");
            GlassProd.writeDataToFile(dsMSDcorP, "MSDcorP.dat");
            GlassProd.writeDataToFile(dsCorVisc, "viscCor.dat");
            for (int i = 0; i < configStorageMSD.getLastConfigIndex() - 1; i++) {
                meterGs.setConfigIndex(i);
                GlassProd.writeDataToFile(meterGs, "Gs_t" + i + ".dat");
                meterGsA.setConfigIndex(i);
                GlassProd.writeDataToFile(meterGsA, "GsA_t" + i + ".dat");
                meterGsB.setConfigIndex(i);
                if (params.nB>0) GlassProd.writeDataToFile(meterGsB, "GsB_t" + i + ".dat");
            }
            for (int i = params.sfacMinInterval; i < sfacMobilityA.length; i++) {
                if (sfacMobilityA[i].acc.getSampleCount() < 2) continue;
                GlassProd.writeDataToFile(sfacMobilityA[i].acc, "sfacMobilityA" + i + ".dat");
                GlassProd.writeDataToFile(sfacMobilityB[i].acc, "sfacMobilityB" + i + ".dat");
                GlassProd.writeDataToFile(sfacMotionX[i].acc, "sfacMotionx" + i + ".dat");
            }
            double fac = sim.box.getBoundary().getBoxSize().getX(0) / (2 * Math.PI);
            for (int j=0; j<2; j++) {
                GlassProd.writeDataToFile(sfcMobilityACor.makeMeter(j), "sfacMobilityACor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfcMobilityBCor.makeMeter(j), "sfacMobilityBCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfcMotionXXCor.makeMeter(j), "sfacMotionXXCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfcMotionXYCor.makeMeter(j), "sfacMotionXYCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfacDensityX.cor.makeMeter(j), "sfacDensityCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfacPackX.cor.makeMeter(j), "sfacPackCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfacABX.cor.makeMeter(j), "sfacABCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfacKEX.cor.makeMeter(j), "sfacKineticCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(sfacStressX.cor.makeMeter(j), "sfacStressCor" + (j+1)+".dat");

                GlassProd.writeDataToFile(dsCorSFacPackMobilityA.makeMeter(j), "sfacPackMobilityACor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacPackMobilityB.makeMeter(j), "sfacPackMobilityBCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacABMobilityA.makeMeter(j), "sfacABMobilityACor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacABMobilityB.makeMeter(j), "sfacABMobilityBCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacMobilityAB.makeMeter(j), "sfacMobilityABCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacKEMobilityA.makeMeter(j), "sfacKEMobilityACor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacKEMobilityB.makeMeter(j), "sfacKEMobilityBCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacPackPackAB.makeMeter(j), "sfacPackPackABCor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacStressMobilityA.makeMeter(j), "sfacStressMobilityACor" + (j+1)+".dat");
                GlassProd.writeDataToFile(dsCorSFacStressMobilityB.makeMeter(j), "sfacStressMobilityBCor" + (j+1)+".dat");
            }

            GlassProd.writeDataToFile(meterCorrelationSelf, "corSelf.dat");
            GlassProd.writeDataToFile(meterCorrelationSelfMagA, "corSelfMagA.dat");
            if (params.nB>0) GlassProd.writeDataToFile(meterCorrelationSelfMagB, "corSelfMagB.dat");
            for (int i = 0; i < correlationSelf2.getNumDt(); i++) {
                CorrelationSelf2.MeterCorrelationSelf2 m = correlationSelf2.makeMeter(i);
                GlassProd.writeDataToFile(m, "corRSelf_t" + i + ".dat");
            }

            for (int i=0; i<meterFs.length; i++) {
                for (int j=0; j<2; j++) GlassProd.writeDataToFile(meterFs[i][j], filenameFs[i][j]);
                GlassProd.writeDataToFile(meterF[i], filenameF[i]);
            }
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc, meterPerc3}, filenamePerc);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(),
                    meterPerc3.makeImmFractionSource()}, filenameImmFrac);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(sim.speciesA.getLeafType()),
                    meterPerc3.makeImmFractionSource(sim.speciesA.getLeafType())}, filenameImmFracA);
            if (params.nB>0) GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(sim.speciesB.getLeafType()),
                    meterPerc3.makeImmFractionSource(sim.speciesB.getLeafType())}, filenameImmFracB);
            GlassProd.writeDataToFile(meterPerc.makePerclationByImmFracSource(), filenameImmFracPerc);
            GlassProd.writeDataToFile(accPerc0, filenamePerc0);
            GlassProd.writeDataToFile(meterQ4, "Q4" + fileTag + ".dat");
            GlassProd.writeDataToFile(meterL, filenameL);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeChi4Source(), meterPerc3.makeChi4Source()}, "chi4Star" + fileTag + ".dat");
            GlassProd.writeDataToFile(meterQ4.makeChi4Meter(), "chi4" + fileTag + ".dat");
            GlassProd.writeDataToFile(meterAlpha2A, filenameAlpha2A);
            GlassProd.writeDataToFile(meterAlpha2B, filenameAlpha2B);
            GlassProd.writeDataToFile(sfacDensity.acc, "sfac.dat");
            GlassProd.writeDataToFile(sfacA.acc, "sfacA.dat");
            GlassProd.writeDataToFile(sfacB.acc, "sfacB.dat");
            GlassProd.writeDataToFile(sfacPack.acc, "sfacPack.dat");
            GlassProd.writeDataToFile(sfacAB.acc, "sfacAB.dat");

            GlassProd.writeDataToFile(meterMSD, meterMSD.getError(), meterMSD.getStdev(), filenameMSD);
            GlassProd.writeDataToFile(dsCorMSD.getFullCorrelation(), "msdFullCor.dat");
            GlassProd.writeDataToFile(dsCorVisc.getFullCorrelation(), "viscFullCor.dat");
            GlassProd.writeDataToFile(meterMSDA, meterMSDA.getError(), meterMSDA.getStdev(), filenameMSDA);
            if (params.nB>0) GlassProd.writeDataToFile(meterMSDB, meterMSDB.getError(), meterMSDB.getStdev(), filenameMSDB);
            GlassProd.writeDataToFile(meterVAC, meterVAC.errData, filenameVAC);
            GlassProd.writeDataToFile(pTensorAccumVisc, pTensorAccumVisc.errData, filenameVisc);
            if (params.doLinear){
                GlassProd.writeDataToFile(meterLinearMSD, meterLinearMSD.errData, filenameLinearMSD);
                GlassProd.writeDataToFile(meterLinearVAC, meterLinearVAC.errData, filenameLinearVAC);
                GlassProd.writeDataToFile(pTensorLinearAccumVisc, pTensorLinearAccumVisc.errData, filenameLinearVisc);
            }
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        long time1 = System.nanoTime();
        double sim_time = (time1 - time0) / 1e9;
        System.out.printf("\ntime: %3.2f s\n", sim_time);
    }

    public static class SimParams extends SimGlass.GlassParams {
        public long numStepsEqIsothermal = 10000;
        public long numStepsIsothermal = 10000;
        public long numSteps = 1000000;
        public double minDrFilter = 0.4;
        public double temperatureMelt = 0;
        public double[] qx = new double[]{7.0}; // for Fs
        public boolean doPxyAutocor = false;
        public int sfacMinInterval = 6;
        public int[] randomSeeds = null;
        public double maxWalltime = Double.POSITIVE_INFINITY;
    }

    public static double[] readMSD(String filename) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            ArrayList<Double> vals = new ArrayList<>();
            String l = null;
            while ((l = br.readLine()) != null) {
                String[] bits = l.split(" ");
                vals.add(Double.parseDouble(bits[1]));
            }
            double[] rv = new double[vals.size()+1];
            for (int i=0; i<vals.size(); i++) {
                rv[i+1] = vals.get(i);
            }
            return rv;
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
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

    public static void writeDataToFile(double[][] data, String filename) throws IOException {
        FileWriter fw = new FileWriter(filename);
        for (int i = 0; i < data.length; i++) {
            for (int j=0; j<data[i].length; j++) {
                double y = data[i][j];
                fw.write(String.format("% 8.5f ", y));
            }
            fw.write("\n");
        }
        fw.close();
    }

    public static void writeCombinedDataToFile(IDataSource[] meters, String filename) {
        List<double[]> allData = new ArrayList<>();
        for (IDataSource meter : meters) {
            IData data = meter.getData();
            IData xData = ((DataFunction.DataInfoFunction) meter.getDataInfo()).getXDataSource().getIndependentData(0);
            for (int i = 0; i < xData.getLength(); i++) {
                double y = data.getValue(i);
                if (Double.isNaN(y)) continue;
                allData.add(new double[]{xData.getValue(i), y});
            }
        }
        allData.sort(new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return Double.compare(o1[0], o2[0]);
            }
        });

        try {
            FileWriter fw = new FileWriter(filename);
            for (double[] xy : allData) {
                fw.write(xy[0] + " " + xy[1] + "\n");
            }
            fw.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static class StepsState implements Statefull {
        public long numStepsEqIsothermal;
        public long numStepsIsothermal;
        public long numSteps;

        @Override
        public void saveState(Writer fw) throws IOException {
            fw.write(numStepsEqIsothermal+" "+numStepsIsothermal+" "+numSteps+"\n");
        }

        @Override
        public void restoreState(BufferedReader br) throws IOException {
            String[] bits = br.readLine().split(" ");
            numStepsEqIsothermal = Long.parseLong(bits[0]);
            numStepsIsothermal = Long.parseLong(bits[1]);
            numSteps = Long.parseLong(bits[2]);
        }
    }
}