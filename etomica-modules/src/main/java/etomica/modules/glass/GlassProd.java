/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

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
            params.numStepsEq = 100000;
            params.numSteps = 1000000;
            params.minDrFilter = 0.4;
            params.qx = 7.0;
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
        double temperature0 = Math.max(params.temperatureMelt, params.temperature);
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
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        tracer.setDataSink(pAccumulator);

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

        //VAC
        DataSourceVAC meterVAC = new DataSourceVAC(configStorageMSD);
        configStorageMSD.addListener(meterVAC);

        //Fs
        DataSourceFs meterFs = new DataSourceFs(configStorageMSD);
        Vector q = sim.getSpace().makeVector();
        q.setX(0, params.qx);
        meterFs.setQ(q);
        configStorageMSD.addListener(meterFs);

        //Percolation
        AtomTestDeviation atomFilterDeviation = new AtomTestDeviation(sim.box, configStorageMSD);
        atomFilterDeviation.setMinDistance(params.minDrFilter);
        atomFilterDeviation.setDoMobileOnly(false);
        DataSourcePercolation meterPerc = new DataSourcePercolation(configStorageMSD, atomFilterDeviation, params.log2StepMin);
        configStorageMSD.addListener(meterPerc);

        AtomTestDeviation atomFilterDeviation3 = new AtomTestDeviation(sim.box, configStorageMSD3);
        atomFilterDeviation3.setMinDistance(params.minDrFilter);
        atomFilterDeviation3.setDoMobileOnly(false);
        DataSourcePercolation meterPerc3 = new DataSourcePercolation(configStorageMSD3, atomFilterDeviation3, params.log2StepMin, meterPerc.getHistogram());
        configStorageMSD3.addListener(meterPerc3);

        DataSourcePercolation0 meterPerc0 = new DataSourcePercolation0(sim.box, sim.getRandom());
        meterPerc0.setImmFracs(new double[]{0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.55, 0.65, 0.75, 0.85, 1});
        AccumulatorAverageFixed accPerc0 = new AccumulatorAverageFixed(10);
        DataPumpListener pumpPerc0 = new DataPumpListener(meterPerc0, accPerc0, 10000);
        sim.integrator.getEventManager().addListener(pumpPerc0);

        DataSourceQ4 meterQ4 = new DataSourceQ4(configStorageMSD, 8);
        meterQ4.setMaxDr(0.2);
        configStorageMSD.addListener(meterQ4);

        //Strings
        DataSourceStrings meterL = new DataSourceStrings(configStorageMSD, 4);
        configStorageMSD.addListener(meterL);

        //Alpha2
        DataSourceAlpha2 meterAlpha2 = new DataSourceAlpha2(configStorageMSD);
        configStorageMSD.addListener(meterAlpha2);

        //S(q)
        AccumulatorAverageFixed accSFac = new AccumulatorAverageFixed(1);  // just average, no uncertainty
        double cut1 = 10;
        if (numAtoms > 500) cut1 /= Math.pow(numAtoms / 500.0, 1.0 / sim.getSpace().D());
        MeterStructureFactor meterSFac = new MeterStructureFactor(sim.box, cut1);
        meterSFac.setNormalizeByN(true);
        DataPumpListener pumpSFac = new DataPumpListener(meterSFac, accSFac, 1000);
        sim.integrator.getEventManager().addListener(pumpSFac);
        double vB = sim.getSpace().powerD(sim.sigmaB);
        ((MeterStructureFactor.AtomSignalSourceByType) meterSFac.getSignalSource()).setAtomTypeFactor(sim.speciesB.getAtomType(0), vB);

        AccumulatorAverageFixed[] accSFacMobility = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMobility signalMobility = new AtomSignalMobility(configStorageMSD);
            signalMobility.setPrevConfig(i);
            MeterStructureFactor meterSFacMobility = new MeterStructureFactor(sim.box, 3, signalMobility);
            meterSFacMobility.setNormalizeByN(true);
            DataFork forkSFacMobility = new DataFork();
            DataPump pumpSFacMobility = new DataPump(meterSFacMobility, forkSFacMobility);
            accSFacMobility[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            forkSFacMobility.addDataSink(accSFacMobility[i]);
            // ensures pump fires when config with delta t is available
            ConfigurationStoragePumper cspMobility = new ConfigurationStoragePumper(pumpSFacMobility, configStorageMSD);
            cspMobility.setPrevStep(Math.max(i, 9));
            configStorageMSD.addListener(cspMobility);
        }

        AccumulatorAverageFixed[] accSFacMotion = new AccumulatorAverageFixed[30];
        for (int i = 0; i < 30; i++) {
            AtomSignalMotion signalMotion = new AtomSignalMotion(configStorageMSD, 0);
            signalMotion.setPrevConfig(i);
            MeterStructureFactor meterSFacMotion = new MeterStructureFactor(sim.box, 3, signalMotion);
            meterSFacMotion.setNormalizeByN(true);
            DataFork forkSFacMotion = new DataFork();
            DataPump pumpSFacMotion = new DataPump(meterSFacMotion, forkSFacMotion);
            accSFacMotion[i] = new AccumulatorAverageFixed(1);  // just average, no uncertainty
            forkSFacMotion.addDataSink(accSFacMotion[i]);
            // ensures pump fires when config with delta t is available
            ConfigurationStoragePumper cspMotion = new ConfigurationStoragePumper(pumpSFacMotion, configStorageMSD);
            cspMotion.setPrevStep(Math.max(i, 9));
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

        AtomSignalStress signalStress0 = new AtomSignalStress(stressSource, 0, 1);
        MeterStructureFactor meterSFacStress0 = new MeterStructureFactor(sim.box, 3, signalStress0);
        meterSFacStress0.setNormalizeByN(true);
        AccumulatorAverageFixed accSFacStress = new AccumulatorAverageFixed(1);
        DataPumpListener pumpSFacStress = new DataPumpListener(meterSFac, accSFacStress, 100);
        if (sim.getSpace().D() == 3) {
            AtomSignalStress signalStress1 = new AtomSignalStress(stressSource, 0, 2);
            MeterStructureFactor meterSFacStress1 = new MeterStructureFactor(sim.box, 3, signalStress1);
            meterSFacStress1.setNormalizeByN(true);
            AtomSignalStress signalStress2 = new AtomSignalStress(stressSource, 1, 2);
            MeterStructureFactor meterSFacStress2 = new MeterStructureFactor(sim.box, 3, signalStress2);
            meterSFacStress2.setNormalizeByN(true);
            MeterStructureFactorStress3 meterStructureFactorStress3 = new MeterStructureFactorStress3(new MeterStructureFactor[]{meterSFacStress0, meterSFacStress1, meterSFacStress2});
            pumpSFacStress = new DataPumpListener(meterStructureFactorStress3, accSFacStress, 100);
        }
        sim.integrator.getEventManager().addListener(pumpSFacStress);

        Vector[] wv = meterSFac.getWaveVectors();
        java.util.List<Vector> myWV = new ArrayList<>();
        double L = sim.box.getBoundary().getBoxSize().getX(0);
        double wvMax2 = 8.01 * Math.PI / L;
        for (Vector vector : wv) {
            int nd = 0;
            for (int i = 0; i < vector.getD(); i++) if (vector.getX(i) != 0) nd++;
            if (vector.squared() > wvMax2 * wvMax2 || nd > 1) continue;
            myWV.add(vector);
        }
        int minIntervalSfac2 = 8;
        wv = myWV.toArray(new Vector[0]);
        MeterStructureFactor[] meterSFacMotion2 = new MeterStructureFactor[30];
        int[] motionMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, 0);
        StructureFactorComponentCorrelation sfcMotionCor = new StructureFactorComponentCorrelation(motionMap, configStorageMSD);
        sfcMotionCor.setMinInterval(minIntervalSfac2);
        MeterStructureFactor[] meterSFacMobility2 = new MeterStructureFactor[30];
        int[] mobilityMap = StructureFactorComponentCorrelation.makeWaveVectorMap(wv, -1);
        StructureFactorComponentCorrelation sfcMobilityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcMobilityCor.setMinInterval(minIntervalSfac2);

        MeterStructureFactor meterSFacDensity2 = new MeterStructureFactor(sim.box, 3);
        meterSFacDensity2.setNormalizeByN(true);
        meterSFacDensity2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcDensityCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcDensityCor.setMinInterval(0);
        DataSinkBlockAveragerSFac dsbaSfacDensity2 = new DataSinkBlockAveragerSFac(configStorageMSD, 0, meterSFacDensity2);
        dsbaSfacDensity2.addSink(sfcDensityCor);
        DataPump pumpSFacDensity2 = new DataPump(meterSFacDensity2, dsbaSfacDensity2);
        ConfigurationStoragePumper cspDensity2 = new ConfigurationStoragePumper(pumpSFacDensity2, configStorageMSD);
        configStorageMSD.addListener(cspDensity2);
        cspDensity2.setPrevStep(0);
        DataSourceCorrelation dsCorSFacDensityMobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacDensity2.addSink(dsCorSFacDensityMobility.makeReceiver(0));

        MeterStructureFactor.AtomSignalSourceByType atomSignalPacking = new MeterStructureFactor.AtomSignalSourceByType();
        atomSignalPacking.setAtomTypeFactor(sim.speciesB.getLeafType(), vB);
        MeterStructureFactor meterSFacPacking2 = new MeterStructureFactor(sim.box, 3, atomSignalPacking);
        meterSFacPacking2.setNormalizeByN(true);
        meterSFacPacking2.setWaveVec(wv);
        StructureFactorComponentCorrelation sfcPackingCor = new StructureFactorComponentCorrelation(mobilityMap, configStorageMSD);
        sfcPackingCor.setMinInterval(0);
        DataSinkBlockAveragerSFac dsbaSfacPacking2 = new DataSinkBlockAveragerSFac(configStorageMSD, 0, meterSFacPacking2);
        dsbaSfacPacking2.addSink(sfcPackingCor);
        DataPump pumpSFacPacking2 = new DataPump(meterSFacPacking2, dsbaSfacPacking2);
        ConfigurationStoragePumper cspPacking2 = new ConfigurationStoragePumper(pumpSFacPacking2, configStorageMSD);
        configStorageMSD.addListener(cspPacking2);
        cspPacking2.setPrevStep(0);
        DataSourceCorrelation dsCorSFacPackingMobility = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacPacking2.addSink(dsCorSFacPackingMobility.makeReceiver(0));
        DataSourceCorrelation dsCorSFacPackingDensity = new DataSourceCorrelation(configStorageMSD, mobilityMap.length);
        dsbaSfacPacking2.addSink(dsCorSFacPackingDensity.makeReceiver(0));
        dsbaSfacDensity2.addSink(dsCorSFacPackingDensity.makeReceiver(1));

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
        }


        double xGsMax = 3;
        int gsMinConfig = 5;
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
        sim.getController().actionPerformed();

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
            double phi;
            if(params.D == 2){
                phi = Math.PI/4*(params.nA+params.nB/(1.4*1.4))/volume;
            }else{
                phi = Math.PI/6*(params.nA+params.nB/(1.4*1.4*1.4))/volume;
            }
            System.out.println("rho: " + params.density + "  phi: " + phi+"\n");
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
            System.out.println("rho: " + params.density+"\n");
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
                if (foo[mobilityMap[j]] != 0) continue;
                foo[mobilityMap[j]] = 1;
                String label = String.format("q%d%d%d.dat", Math.round(Math.abs(myWV.get(j).getX(0)) * fac),
                        Math.round(Math.abs(myWV.get(j).getX(1)) * fac),
                        Math.round(Math.abs(myWV.get(j).getX(2)) * fac));

                StructureFactorComponentCorrelation.Meter m = sfcMobilityCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacMobilityCor_" + label);
                m = sfcDensityCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacDensityCor_" + label);
                m = sfcPackingCor.makeMeter(mobilityMap[j]);
                GlassProd.writeDataToFile(m, "sfacPackingCor_" + label);
                DataSourceCorrelation.Meter mm = dsCorSFacDensityMobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacDensityMobilityCor_" + label);
                mm = dsCorSFacPackingMobility.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacPackingMobilityCor_" + label);
                mm = dsCorSFacPackingDensity.makeMeter(j);
                GlassProd.writeDataToFile(mm, "sfacPackingDensityCor_" + label);
            }

            foo = new int[motionMap.length];
            for (int j = 0; j < motionMap.length; j++) {
                if (foo[motionMap[j]] != 0) continue;
                foo[motionMap[j]] = 1;
                String label = String.format("q%d%d%d.dat", Math.round(myWV.get(j).getX(0) * fac),
                        Math.round(myWV.get(j).getX(1) * fac),
                        Math.round(myWV.get(j).getX(2) * fac));

                StructureFactorComponentCorrelation.Meter m = sfcMotionCor.makeMeter(motionMap[j]);
                GlassProd.writeDataToFile(m, "sfacMotionCor_" + label);
            }

            GlassProd.writeDataToFile(accSFacStress, "sfacStress.dat");
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
            GlassProd.writeDataToFile(meterPerc0, filenamePerc0);
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
        public int log2StepMin = 9;
        public double temperatureMelt = 0;
        public double qx = 7.0;
        public boolean doPxyAutocor = false;
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
