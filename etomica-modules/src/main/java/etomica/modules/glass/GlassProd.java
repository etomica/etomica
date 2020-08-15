/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.activity.ActivityIntegrate;
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

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

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
            params.numStepsEq = 10000;
            params.numSteps = 100000;
            params.minDrFilter = 0.4;
            params.qx = 7.0;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep, params.randomSeeds);
        if (params.randomSeeds == null) System.out.println("random seeds: " + Arrays.toString(sim.getRandomSeeds()));
        else System.out.println("set random seeds: " + Arrays.toString(params.randomSeeds));
        System.out.println(params.D + "D " + sim.potentialChoice);
        System.out.println("nA:nB = " + params.nA + ":" + params.nB);
        int numAtoms = params.nA + params.nB;
        double rho = params.density;
        System.out.println("T = " + params.temperature);
        System.out.println(params.numSteps + " MD steps after " + params.numStepsEq + " equilibaration steps , using dt = " + sim.integrator.getTimeStep());
        double volume = sim.box.getBoundary().volume();
        if (params.potential == SimGlass.PotentialChoice.HS) {
            double phi;
            if (params.D == 2) {
                phi = Math.PI / 4 * (params.nA + params.nB / (1.4 * 1.4)) / volume;
            } else {
                phi = Math.PI / 6 * (params.nA + params.nB / (1.4 * 1.4 * 1.4)) / volume;
            }
            System.out.println("rho: " + params.density + "  phi: " + phi + "\n");
        } else {
            System.out.println("rho: " + params.density + "\n");
        }
        sim.initConfig();


        //Equilibration
        double temperature0 = Math.max(params.temperatureMelt, params.temperature);
        if (temperature0 > params.temperature) System.out.println("Equilibrating at T=" + temperature0);
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(temperature0);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numStepsEq / 2);

        AccumulatorAverageFixed accE = new AccumulatorAverageFixed(1);
        MeterEnergy meterE = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        DataPumpListener pumpE = new DataPumpListener(meterE, accE, 10);
        if (sim.potentialChoice != SimGlass.PotentialChoice.HS && temperature0 == params.temperature) {
            sim.integrator.getEventManager().addListener(pumpE);
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numStepsEq / 2);

        if (temperature0 > params.temperature) {
            System.out.println("Equilibrating at T=" + params.temperature);
            sim.integrator.setTemperature(params.temperature);
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numStepsEq / 2);

            if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
                sim.integrator.getEventManager().addListener(pumpE);
            }

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numStepsEq / 2);
        }

        if (sim.potentialChoice != SimGlass.PotentialChoice.HS) {
            // Try to use the average total energy from the 2nd half of equilibration
            // to set the current energy now (via the temperature) so that the average
            // temperature during production will be approximately equal to the set
            // temperature.
            sim.integrator.getEventManager().removeListener(pumpE);
            double avgE = accE.getData(accE.AVERAGE).getValue(0);
            System.out.println("average energy during second half of eq: " + avgE / numAtoms);
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotentialMaster(), sim.box);
            double pe = meterPE.getDataAsScalar();
//            double ke = new MeterKineticEnergy(sim.box).getDataAsScalar();
//            System.out.println("kinetic energy at end of eq: "+ke/numAtoms);
//            System.out.println("potential energy at end of eq: "+pe/numAtoms);
            double newKE = avgE - pe;
            double nowTemp = newKE * 2.0 / params.D / numAtoms;
            double oldTemp = new MeterTemperature(sim.box, params.D).getDataAsScalar();
            System.out.println("setting temp " + oldTemp + " => " + nowTemp); // + " (ke "+ke+" => " + newKE / numAtoms + ")");
            sim.integrator.setIntegratorMC(null, 0);
            sim.integrator.setIsothermal(false);
            sim.integrator.setTemperature(nowTemp);
//            double newNewKE = new MeterKineticEnergy(sim.box).getDataAsScalar();
//            System.out.println("actual KE => "+newNewKE/numAtoms);
        }

        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);

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


        AccumulatorAverageFixed gTensorAccumulator = new AccumulatorAverageFixed(blocksize);
        if (params.potential != SimGlass.PotentialChoice.HS) {
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
        boolean doPxyAutocor = params.doPxyAutocor && sim.potentialChoice != SimGlass.PotentialChoice.HS;
        if (doPxyAutocor) {
            pTensorFork.addDataSink(dpxyAutocor);
        }
        dpxyAutocor.setPushInterval(16384);


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
        configStorageMSD.addListener(meterMSDB);

        DataSourceCorMSD dsCorMSD = new DataSourceCorMSD(sim.integrator);
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

        //VAC
        configStorageMSD.setDoVelocity(true);
        DataSourceVAC meterVAC = new DataSourceVAC(configStorageMSD);
        configStorageMSD.addListener(meterVAC);

        //Fs
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        Vector q = sim.getSpace().makeVector();
        q.setX(0, params.qx);
        meterFs.setQ(q);
        configStorageMSD.addListener(meterFs);

        //Percolation
        int percMinInterval = 11;
        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorageMSD);
        atomFilterDeviation.setMinDistance(params.minDrFilter);
        atomFilterDeviation.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSD, atomFilterDeviation, percMinInterval);
        configStorageMSD.addListener(meterPerc);

        AtomTestDeviation atomFilterDeviation3 = new AtomTestDeviation(sim.box, configStorageMSD3);
        atomFilterDeviation3.setMinDistance(params.minDrFilter);
        atomFilterDeviation3.setDoMobileOnly(false);
        DataSourcePercolation meterPerc3 = new DataSourcePercolation(configStorageMSD3, atomFilterDeviation3, percMinInterval - 1, meterPerc.getHistogram());
        configStorageMSD3.addListener(meterPerc3);

        DataSourcePercolation0 meterPerc0 = new DataSourcePercolation0(sim.box, sim.getRandom());
        meterPerc0.setImmFracs(new double[]{0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27});
        AccumulatorAverageFixed accPerc0 = new AccumulatorAverageFixed(10);
        DataPumpListener pumpPerc0 = new DataPumpListener(meterPerc0, accPerc0, 10000);
        sim.integrator.getEventManager().addListener(pumpPerc0);

        DataSourceQ4 meterQ4 = new DataSourceQ4(configStorageMSD, percMinInterval + 1);
        meterQ4.setMaxDr(0.2);
        configStorageMSD.addListener(meterQ4);

        //Strings
        DataSourceStrings meterL = new DataSourceStrings(configStorageMSD, 6);
        configStorageMSD.addListener(meterL);

        //Alpha2
        DataSourceAlpha2 meterAlpha2 = new DataSourceAlpha2(configStorageMSD);
        configStorageMSD.addListener(meterAlpha2);

        //S(q)
        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);  // just average, no uncertainty
        double cut10 = 10;
        if (numAtoms > 500) cut10 /= Math.pow(numAtoms / 500.0, 1.0 / sim.getSpace().D());
        MeterStructureFactor meterSFac = new MeterStructureFactor(sim.box, cut10);
        meterSFac.setNormalizeByN(true);
        DataPumpListener pumpSFac = new DataPumpListener(meterSFac, accSFac, 5000);
        sim.integrator.getEventManager().addListener(pumpSFac);
        double vB = sim.getSpace().powerD(sim.sigmaB);
        ((MeterStructureFactor.AtomSignalSourceByType) meterSFac.getSignalSource()).setAtomTypeFactor(sim.speciesB.getAtomType(0), vB);

        double cut3 = 3;
        if (numAtoms > 500) cut3 /= Math.pow(numAtoms / 500.0, 1.0 / sim.getSpace().D());
        AccumulatorAverageFixed[] accSFacMobility = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMobility signalMobility = new AtomSignalMobility(configStorageMSD);
            signalMobility.setPrevConfig(i + 1);
            MeterStructureFactor meterSFacMobility = new MeterStructureFactor(sim.box, cut3, signalMobility);
            meterSFacMobility.setNormalizeByN(true);
            accSFacMobility[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            DataPump pumpSFacMobility = new DataPump(meterSFacMobility, accSFacMobility[i]);
            // ensures pump fires when config with delta t is available
            ConfigurationStoragePumper cspMobility = new ConfigurationStoragePumper(pumpSFacMobility, configStorageMSD);
            cspMobility.setPrevStep(Math.max(i, 11));
            configStorageMSD.addListener(cspMobility);
        }

        AccumulatorAverageFixed[] accSFacMotion = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMotion signalMotion = new AtomSignalMotion(configStorageMSD, 0);
            signalMotion.setPrevConfig(i + 1);
            MeterStructureFactor meterSFacMotion = new MeterStructureFactor(sim.box, cut3, signalMotion);
            meterSFacMotion.setNormalizeByN(true);
            accSFacMotion[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            DataPump pumpSFacMotion = new DataPump(meterSFacMotion, accSFacMotion[i]);
            // ensures pump fires when config with delta t is available
            ConfigurationStoragePumper cspMotion = new ConfigurationStoragePumper(pumpSFacMotion, configStorageMSD);
            cspMotion.setPrevStep(Math.max(i, 11));
            configStorageMSD.addListener(cspMotion);
        }

        AtomStressSource stressSource = null;
        if (sim.potentialChoice == SimGlass.PotentialChoice.HS) {
            AtomHardStressCollector ahsc = new AtomHardStressCollector((IntegratorHard) sim.integrator);
            ((IntegratorHard) sim.integrator).addCollisionListener(ahsc);
            stressSource = ahsc;
        } else {
            PotentialCalculationForceSumGlass pcForce = new PotentialCalculationForceSumGlass(sim.box);
            ((IntegratorVelocityVerlet) sim.integrator).setForceSum(pcForce);
            stressSource = pcForce;
        }
        int[][] normalComps = new int[sim.getSpace().getD()][2];
        for (int i = 0; i < normalComps.length; i++) {
            normalComps[i][0] = normalComps[i][1] = i;
        }
        AtomSignalStress signalStressNormal = new AtomSignalStress(stressSource, normalComps);

        Vector[] wv = meterSFac.getWaveVectors();
        List<Vector> myWV = new ArrayList<>();
        double L = sim.box.getBoundary().getBoxSize().getX(0);
        double wvMax2 = 4.01 * Math.PI / L;
        for (Vector vector : wv) {
            int nd = 0;
            for (int i = 0; i < vector.getD(); i++) if (vector.getX(i) != 0) nd++;
            if (vector.squared() > wvMax2 * wvMax2 || nd > 1) continue;
            myWV.add(vector);
        }
        int minIntervalSfac2 = 6;
        wv = myWV.toArray(new Vector[0]);
        int[] motionMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, 0);
        int[] mobilityMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, -1);

        MeterStructureFactor meterSFacStress2 = new MeterStructureFactor(sim.box, 3, signalStressNormal);
        meterSFacStress2.setNormalizeByN(true);
        meterSFacStress2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcStress2Cor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcStress2Cor.setMinInterval(params.sfacMinInterval);
        DataSinkBlockAveragerSFac dsbaSfacStress2 = new DataSinkBlockAveragerSFac(configStorageMSD, params.sfacMinInterval, meterSFacStress2);
        dsbaSfacStress2.addSink(sfcStress2Cor);
        DataPump pumpSFacStress2Cor = new DataPump(meterSFacStress2, dsbaSfacStress2);
        ConfigurationStoragePumper cspStress2 = new ConfigurationStoragePumper(pumpSFacStress2Cor, configStorageMSD);
        configStorageMSD.addListener(cspStress2);
        cspStress2.setPrevStep(params.sfacMinInterval);
        DataSourceCorrelation dsCorSFacStress2Mobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacStress2.addSink(dsCorSFacStress2Mobility.makeReceiver(0));

        MeterStructureFactor[] meterSFacMotion2 = new MeterStructureFactor[30];
        StructureFactorComponentCorrelation sfcMotionCor = new StructureFactorComponentCorrelation(motionMap, configStorageMSD);
        sfcMotionCor.setMinInterval(params.sfacMinInterval);
        MeterStructureFactor[] meterSFacMobility2 = new MeterStructureFactor[30];
        StructureFactorComponentCorrelation sfcMobilityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcMobilityCor.setMinInterval(params.sfacMinInterval);

        MeterStructureFactor meterSFacDensity2 = new MeterStructureFactor(sim.box, 3);
        meterSFacDensity2.setNormalizeByN(true);
        meterSFacDensity2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcDensityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcDensityCor.setMinInterval(params.sfacMinInterval);
        DataSinkBlockAveragerSFac dsbaSfacDensity2 = new DataSinkBlockAveragerSFac(configStorageMSD, params.sfacMinInterval, meterSFacDensity2);
        dsbaSfacDensity2.addSink(sfcDensityCor);
        DataPump pumpSFacDensity2 = new DataPump(meterSFacDensity2, dsbaSfacDensity2);
        ConfigurationStoragePumper cspDensity2 = new ConfigurationStoragePumper(pumpSFacDensity2, configStorageMSD);
        configStorageMSD.addListener(cspDensity2);
        cspDensity2.setPrevStep(params.sfacMinInterval);
        DataSourceCorrelation dsCorSFacDensityMobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacDensity2.addSink(dsCorSFacDensityMobility.makeReceiver(0));

        MeterStructureFactor.AtomSignalSourceByType atomSignalPacking = new MeterStructureFactor.AtomSignalSourceByType();
        atomSignalPacking.setAtomTypeFactor(sim.speciesB.getLeafType(), sim.potentialChoice == SimGlass.PotentialChoice.HS ? vB : 0);
        MeterStructureFactor meterSFacPacking2 = new MeterStructureFactor(sim.box, 3, atomSignalPacking);
        meterSFacPacking2.setNormalizeByN(true);
        meterSFacPacking2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcPackingCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcPackingCor.setMinInterval(params.sfacMinInterval);
        DataSinkBlockAveragerSFac dsbaSfacPacking2 = new DataSinkBlockAveragerSFac(configStorageMSD, params.sfacMinInterval, meterSFacPacking2);
        dsbaSfacPacking2.addSink(sfcPackingCor);
        DataPump pumpSFacPacking2 = new DataPump(meterSFacPacking2, dsbaSfacPacking2);
        ConfigurationStoragePumper cspPacking2 = new ConfigurationStoragePumper(pumpSFacPacking2, configStorageMSD);
        configStorageMSD.addListener(cspPacking2);
        cspPacking2.setPrevStep(params.sfacMinInterval);
        DataSourceCorrelation dsCorSFacPackingMobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacPacking2.addSink(dsCorSFacPackingMobility.makeReceiver(0));
        DataSourceCorrelation dsCorSFacPackingDensity = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacPacking2.addSink(dsCorSFacPackingDensity.makeReceiver(0));
        dsbaSfacDensity2.addSink(dsCorSFacPackingDensity.makeReceiver(1));

        AtomSignalKineticEnergy atomSignalKE = new AtomSignalKineticEnergy();
        atomSignalKE.setDoSubtractAvg(1.5 * sim.integrator.getTemperature());
        MeterStructureFactor meterSFacKE = new MeterStructureFactor(sim.box, 3, atomSignalKE);
        meterSFacKE.setNormalizeByN(true);
        meterSFacKE.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcKECor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcKECor.setMinInterval(params.sfacMinInterval);
        DataSinkBlockAveragerSFac dsbaSfacKE = new DataSinkBlockAveragerSFac(configStorageMSD, params.sfacMinInterval, meterSFacKE);
        dsbaSfacKE.addSink(sfcKECor);
        DataPump pumpSFacKECor = new DataPump(meterSFacKE, dsbaSfacKE);
        ConfigurationStoragePumper cspKE = new ConfigurationStoragePumper(pumpSFacKECor, configStorageMSD);
        configStorageMSD.addListener(cspKE);
        cspKE.setPrevStep(params.sfacMinInterval);
        DataSourceCorrelation dsCorSFacKEMobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacKE.addSink(dsCorSFacKEMobility.makeReceiver(0));

        for (int i = 0; i < 30; i++) {
            AtomSignalMotion signalMotion = new AtomSignalMotion(configStorageMSD, 0);
            signalMotion.setPrevConfig(i + 1);
            meterSFacMotion2[i] = new MeterStructureFactor(sim.box, 3, signalMotion);
            meterSFacMotion2[i].setNormalizeByN(true);
            meterSFacMotion2[i].setWaveVec(wv);
            DataPump pumpSFacMotion2 = new DataPump(meterSFacMotion2[i], sfcMotionCor.makeSink(i, meterSFacMotion2[i]));
            ConfigurationStoragePumper cspMotion2 = new ConfigurationStoragePumper(pumpSFacMotion2, configStorageMSD);
            cspMotion2.setPrevStep(i);
            cspMotion2.setBigStep(minIntervalSfac2);
            configStorageMSD.addListener(cspMotion2);

            AtomSignalMobility signalMobility = new AtomSignalMobility(configStorageMSD);
            signalMobility.setPrevConfig(i + 1);
            meterSFacMobility2[i] = new MeterStructureFactor(sim.box, 3, signalMobility);
            meterSFacMobility2[i].setNormalizeByN(true);
            meterSFacMobility2[i].setWaveVec(wv);
            DataFork sfacMobility2Fork = new DataFork();
            sfacMobility2Fork.addDataSink(sfcMobilityCor.makeSink(i, meterSFacMobility2[i]));
            DataPump pumpSFacMobility2 = new DataPump(meterSFacMobility2[i], sfacMobility2Fork);
            ConfigurationStoragePumper cspMobility2 = new ConfigurationStoragePumper(pumpSFacMobility2, configStorageMSD);
            cspMobility2.setPrevStep(i);
            cspMobility2.setBigStep(minIntervalSfac2);
            configStorageMSD.addListener(cspMobility2);
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacDensityMobility));
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacPackingMobility));
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacStress2Mobility));
            sfacMobility2Fork.addDataSink(new StructorFactorComponentExtractor(meterSFacMobility2, i, dsCorSFacKEMobility));
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
        configStorageMSD.addListener(meterGsB);
        meterGsB.getXDataSource().setXMax(xGsMax);

        MeterCorrelationSelf meterCorrelationSelf = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.TOTAL);
        configStorageMSD.addListener(meterCorrelationSelf);
        MeterCorrelationSelf meterCorrelationSelfMagA = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagA.setAtomType(sim.speciesA.getLeafType());
        configStorageMSD.addListener(meterCorrelationSelfMagA);
        MeterCorrelationSelf meterCorrelationSelfMagB = new MeterCorrelationSelf(configStorageMSD, MeterCorrelationSelf.CorrelationType.MAGNITUDE);
        meterCorrelationSelfMagB.setAtomType(sim.speciesB.getLeafType());
        configStorageMSD.addListener(meterCorrelationSelfMagB);

        CorrelationSelf2 correlationSelf2 = new CorrelationSelf2(configStorageMSD, CorrelationSelf2.CorrelationType.TOTAL, 0.001, 20);
        configStorageMSD.addListener(correlationSelf2);

        sim.integrator.getEventManager().addListener(configStorageMSD3);
        sim.integrator.getEventManager().addListener(configStorageMSD);

        //Run
        long time0 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps);

        //Pressure
        DataGroup dataP = (DataGroup)pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg  = dataPAvg.getValue(0);
        double pErr  = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);

        String filenameVisc, filenameMSD, filenameMSDA, filenameMSDB, filenameFs, filenamePerc,
                filenamePerc0, filenameImmFrac, filenameImmFracA, filenameImmFracB, filenameImmFracPerc, filenameL, filenameAlpha2, filenameSq, filenameVAC;

        String fileTag = "";
        String filenamePxyAC = "";
        if(sim.potentialChoice == SimGlass.PotentialChoice.HS){
            System.out.println("Z: " + pAvg / params.density / params.temperature + "  " + pErr / params.density / params.temperature + "  cor: " + pCorr);
            fileTag = String.format("Rho%1.3f", rho);
            filenameVisc = String.format("viscRho%1.3f.out", rho);
            filenameMSD = String.format("msdRho%1.3f.out", rho);
            filenameMSDA = String.format("msdARho%1.3f.out", rho);
            filenameMSDB = String.format("msdBRho%1.3f.out", rho);
            filenameVAC = String.format("vacRho%1.3f.out", rho);
            filenameFs = String.format("fsRho%1.3fQ%1.2f.out", rho, params.qx);
            filenamePerc = String.format("percRho%1.3f.out", rho);
            filenamePerc0 = String.format("perc0Rho%1.3f.out", rho);
            filenameL = String.format("lRho%1.3f.out", rho);
            filenameImmFrac = String.format("immFracRho%1.3f.out", rho);
            filenameImmFracA = String.format("immFracARho%1.3f.out", rho);
            filenameImmFracB = String.format("immFracBRho%1.3f.out", rho);
            filenameImmFracPerc = String.format("immFracPercRho%1.3f.out", rho);
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
            System.out.println("T: " + tAvg +"  "+ tErr +"  cor: "+tCorr);
            System.out.println("Z: " + pAvg/params.density/tAvg +"  "+ pErr/params.density/tAvg  +"  cor: "+pCorr);
            System.out.println("U: " + uAvg / numAtoms + "  " + uErr / numAtoms + "  cor: " + uCorr);
            fileTag = String.format("Rho%1.3fT%1.3f", rho, params.temperature);
            filenameVisc = String.format("viscRho%1.3fT%1.3f.out", rho, params.temperature);
            filenameMSD = String.format("msdRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameMSDA = String.format("msdARho%1.3fT%1.3f.out", rho, params.temperature);
            filenameMSDB = String.format("msdBRho%1.3fT%1.3f.out", rho, params.temperature);
            filenameVAC = String.format("vacRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameFs = String.format("fsRho%1.3fT%1.3fQ%1.2f.out", rho, params.temperature, params.qx);
            filenamePerc = String.format("percRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePerc0 = String.format("perc0Rho%1.3fT%1.3f.out", rho, params.temperature);
            filenameL = String.format("lRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameImmFrac = String.format("immFracRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameImmFracA = String.format("immFracARho%1.3fT%1.3f.out", rho, params.temperature);
            filenameImmFracB = String.format("immFracBRho%1.3fT%1.3f.out", rho, params.temperature);
            filenameImmFracPerc = String.format("immFracPercRho%1.3fT%1.3f.out", rho, params.temperature);
            filenameAlpha2 = String.format("alpha2Rho%1.3fT%1.3f.out",  rho, params.temperature);
            filenamePxyAC = String.format("acPxyRho%1.3fT%1.3f.out",  rho, params.temperature);
            filenameSq = String.format("sqRho%1.3fT%1.3f.out",  rho, params.temperature);
            double V = sim.box.getBoundary().volume();
            System.out.println("G: " + V * avgG / tAvg + " " + V * errG / tAvg + " cor: " + corG + "\n");
        }
        System.out.println("P: " + pAvg +"  "+ pErr +"  cor: "+pCorr);

        try {
            if (doPxyAutocor) {
                GlassProd.writeDataToFile(dpxyAutocor, filenamePxyAC);
            }
            GlassProd.writeDataToFile(dsCorMSD, "MSDcor.dat");
            GlassProd.writeDataToFile(dsCorP, "Pcor.dat");
            GlassProd.writeDataToFile(dsMSDcorP, "MSDcorP.dat");
            for (int i = 0; i < configStorageMSD.getLastConfigIndex() - 1; i++) {
                meterGs.setConfigIndex(i);
                GlassProd.writeDataToFile(meterGs, "Gs_t" + i + ".dat");
                meterGsA.setConfigIndex(i);
                GlassProd.writeDataToFile(meterGsA, "GsA_t" + i + ".dat");
                meterGsB.setConfigIndex(i);
                GlassProd.writeDataToFile(meterGsB, "GsB_t" + i + ".dat");
            }
            for (int i = 0; i < accSFacMobility.length; i++) {
                if (accSFacMobility[i].getSampleCount() < 2) continue;
                GlassProd.writeDataToFile(accSFacMobility[i], "sfacMobility" + i + ".dat");
                GlassProd.writeDataToFile(accSFacMotion[i], "sfacMotionx" + i + ".dat");
            }
            double fac = L / (2 * Math.PI);
            int[] foo = new int[mobilityMap.length];
            for (int j = 0; j < mobilityMap.length; j++) {
                String packing = sim.potentialChoice == SimGlass.PotentialChoice.HS ? "Packing" : "DensityA";
                if (foo[mobilityMap[j]] != 0) continue;
                foo[mobilityMap[j]] = 1;
                String label = String.format("q%d%d%d.dat", Math.round(Math.abs(wv[j].getX(0)) * fac),
                        Math.round(Math.abs(wv[j].getX(1)) * fac),
                        Math.round(Math.abs(wv[j].getX(2)) * fac));

                StructureFactorComponentCorrelation.Meter m = sfcMobilityCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacMobilityCor_" + label);
                m = sfcDensityCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacDensityCor_" + label);
                m = sfcPackingCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfac" + packing + "Cor_" + label);
                m = sfcStress2Cor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacStressCor_" + label);
                m = sfcKECor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacKineticCor_" + label);
                DataSourceCorrelation.Meter mm = dsCorSFacDensityMobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacDensityMobilityCor_" + label);
                mm = dsCorSFacPackingMobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfac" + packing + "MobilityCor_" + label);
                mm = dsCorSFacPackingDensity.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfac" + packing + "DensityCor_" + label);
                mm = dsCorSFacStress2Mobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacStressMobilityCor_" + label);
                mm = dsCorSFacKEMobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacKineticMobilityCor_" + label);
            }

            foo = new int[motionMap.length];
            for (int j = 0; j < motionMap.length; j++) {
                if (foo[motionMap[j]] != 0) continue;
                foo[motionMap[j]] = 1;
                String label = String.format("q%d%d%d.dat", Math.round(wv[j].getX(0) * fac),
                        Math.round(wv[j].getX(1) * fac),
                        Math.round(wv[j].getX(2) * fac));

                StructureFactorComponentCorrelation.Meter m = sfcMotionCor.makeMeter(motionMap[j]);
                GlassProd.writeDataToFile(m, "sfacMotionCor_" + label);
            }

            GlassProd.writeDataToFile(meterCorrelationSelf, "corSelf.dat");
            GlassProd.writeDataToFile(meterCorrelationSelfMagA, "corSelfMagA.dat");
            GlassProd.writeDataToFile(meterCorrelationSelfMagB, "corSelfMagB.dat");
            for (int i = 0; i < correlationSelf2.getNumDt(); i++) {
                CorrelationSelf2.MeterCorrelationSelf2 m = correlationSelf2.makeMeter(i);
                GlassProd.writeDataToFile(m, "corRSelf_t" + i + ".dat");
            }

            GlassProd.writeDataToFile(meterFs, filenameFs);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc, meterPerc3}, filenamePerc);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(),
                    meterPerc3.makeImmFractionSource()}, filenameImmFrac);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(sim.speciesA.getLeafType()),
                    meterPerc3.makeImmFractionSource(sim.speciesA.getLeafType())}, filenameImmFracA);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeImmFractionSource(sim.speciesB.getLeafType()),
                    meterPerc3.makeImmFractionSource(sim.speciesB.getLeafType())}, filenameImmFracB);
            GlassProd.writeDataToFile(meterPerc.makePerclationByImmFracSource(), filenameImmFracPerc);
            GlassProd.writeDataToFile(accPerc0, filenamePerc0);
            GlassProd.writeDataToFile(meterQ4, "Q4" + fileTag + ".out");
            GlassProd.writeDataToFile(meterL, filenameL);
            GlassProd.writeCombinedDataToFile(new IDataSource[]{meterPerc.makeChi4Source(), meterPerc3.makeChi4Source()}, "chi4Star" + fileTag + ".out");
            GlassProd.writeDataToFile(meterQ4.makeChi4Meter(), "chi4" + fileTag + ".out");
            GlassProd.writeDataToFile(meterAlpha2, filenameAlpha2);
            GlassProd.writeDataToFile(meterSFac, filenameSq);

            GlassProd.writeDataToFile(meterMSD, meterMSD.errData, filenameMSD);
            GlassProd.writeDataToFile(meterMSDA, meterMSDA.errData, filenameMSDA);
            GlassProd.writeDataToFile(meterMSDB, meterMSDB.errData, filenameMSDB);
            GlassProd.writeDataToFile(meterVAC, meterVAC.errData, filenameVAC);
            GlassProd.writeDataToFile(pTensorAccumVisc, pTensorAccumVisc.errData, filenameVisc);
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
        long time1 = System.nanoTime();
        double sim_time = (time1 - time0) / 1e9;
        System.out.println(String.format("\ntime: %3.2f s", sim_time));
    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;
        public double minDrFilter = 0.4;
        public double temperatureMelt = 0;
        public double qx = 7.0;
        public boolean doPxyAutocor = false;
        public int sfacMinInterval = 6;
        public int[] randomSeeds = null;
    }

    public static void writeDataToFile(IDataSource meter, String filename) throws IOException {
        writeDataToFile(meter, null, filename);
    }

    public static void writeDataToFile(IDataSource meter, IData errData, String filename) throws IOException {
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
            if (errData == null) {
                fw.write(xData.getValue(i) + " " + y + "\n");
            } else {
                fw.write(xData.getValue(i) + " " + y + " " + errData.getValue(i) + "\n");
            }
        }
        fw.close();
    }

    public static void writeCombinedDataToFile(IDataSource[] meters, String filename) {
        List<double[]> allData = new ArrayList<>();
        for (int j = 0; j < meters.length; j++) {
            IData data = meters[j].getData();
            IData xData = ((DataFunction.DataInfoFunction) meters[j].getDataInfo()).getXDataSource().getIndependentData(0);
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

}
