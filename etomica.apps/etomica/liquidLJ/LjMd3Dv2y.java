/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import java.io.File;

import etomica.action.WriteConfigurationBinary;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IFunction;
import etomica.api.IIntegrator;
import etomica.config.ConfigurationFileBinary;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsingLog;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.data.DataFork;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceScalar;
import etomica.data.DataSplitter;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential0Lrc;
import etomica.space3d.Space3D;
import etomica.units.Energy;
import etomica.units.Null;
import etomica.units.SimpleUnit;
import etomica.util.Function;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.HistoryCollapsingDiscard;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3Dv2y {
    
    public static void main(String[] args) {

        // according to http://link.aip.org/link/doi/10.1063/1.2753149
        // triple point
        // T = 0.694
        // liquid density = 0.845435
        
        // Agrawal and Kofke:
        //      small      large
        // T    0.698    0.687(4)
        // p    0.0013   0.0011
        // rho  0.854    0.850
        
        // http://link.aip.org/link/doi/10.1063/1.4758698 says
        // T = 0.7085(5)
        // P = 0.002264(17)
        // rhoL = 0.8405(3)
        // rhoFCC = 0.9587(2)
        // rhoV = 0.002298(18)

        LjMd3DParams params = new LjMd3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 500;
            params.steps = 100000;
            params.v2 = 0.8;
            params.rcShort = 2.5*Math.pow(params.v2, 1.0/6.0);
            params.y = 1.5;
            params.hybridInterval = 100;
        }

        final int numAtoms = params.numAtoms;
        final double y = params.y;
        final double v2 = params.v2;
        final double density = 1.0/Math.sqrt(v2);
        final double temperature = 4/(y*v2*v2);
        long steps = params.steps;
        double tStep = 0.002*Math.sqrt(y)*v2 * params.tStepFac;
        double rcShort = params.rcShort;
        boolean graphics = params.graphics;
        final int hybridInterval = params.hybridInterval;
        int nAccBlocks = params.nAccBlocks;
        boolean calcMu = params.calcMu;

    	int fastInterval = hybridInterval;
        int longInterval = (numAtoms / 200 / hybridInterval) * hybridInterval;
        if (longInterval == 0) longInterval = hybridInterval;

        int numFastInsert = (numAtoms*hybridInterval)/10;
        if (numFastInsert == 0) numFastInsert = 1;
        int numFullInsert = (numAtoms*hybridInterval)/1000;
        if (numFullInsert < 10) numFastInsert = 10;

    	if (!graphics) {
            System.out.println("Running LJ MD with N="+numAtoms+" at y="+y+" v2="+v2);
    	    System.out.println("  T="+temperature+" density="+density);
    	    System.out.println("time step: "+tStep);
    	    System.out.println(steps+" steps ("+(steps*tStep)+" time units)");
    	    System.out.println("short cutoff: "+rcShort);
    	    System.out.println("hybrid MC interval: "+hybridInterval);

    	    if (calcMu) {
    	        System.out.println("fast Widom insertions: "+numFastInsert+" / 1"); //   full Widom insertions: "+numFullInsert+" / "+longInterval);
    	    }
            else {
                System.out.println("fast energy interval: "+fastInterval+"  full energy interval: "+longInterval);
    	    }
    	}

    	double L = Math.pow(numAtoms/density, 1.0/3.0);
        final LjMd3D sim = new LjMd3D(numAtoms, temperature, density, Double.NaN, tStep, rcShort, 0.494*L, hybridInterval, null);
        
    	
    	if (Double.parseDouble(String.format("%4.2f", v2)) != v2) {
    	    throw new RuntimeException(String.format("you're just trying to cause trouble, use v2=%4.2f", v2));
    	}
    	String configFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms, v2, y);
        File inConfigFile = new File(configFilename+".pos");
        if (inConfigFile.exists()) {
            ConfigurationFileBinary configFile = new ConfigurationFileBinary(configFilename);
            configFile.initializeCoordinates(sim.box);
            System.out.println("Continuing from previous configuration");
        }
        else {
            // try 8x fewer atoms
            String tmpConfigFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms/8, v2, y);
            inConfigFile = new File(tmpConfigFilename+".pos");
            if (inConfigFile.exists()) {
                System.out.println("bringing configuration up from N="+numAtoms/8);
                ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
                ConfigurationFileBinary.replicate(configFile, sim.box, new int[]{2,2,2}, Space3D.getInstance());
            }
            else {
                // try lower temperature, then higher temperature
                boolean success = false;
                for (int i=10; i>=-10; i--) {
                    if (i==0) continue;
                    tmpConfigFilename = String.format("configN%d_V%4.2f_Y%4.2f", numAtoms, v2, y+0.01*i);
                    inConfigFile = new File(tmpConfigFilename+".pos");
                    if (inConfigFile.exists()) {
                        System.out.println("bringing configuration from y="+(y+0.01*i));
                        success = true;
                        ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
                        configFile.initializeCoordinates(sim.box);
                        break;
                    }
                }
                if (!success) {
                    // try a different density (higher density first)
                    for (int i=-10; i>=10; i--) {
                        if (i==0) continue;
                        tmpConfigFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms, v2+0.01*i, y);
                        inConfigFile = new File(tmpConfigFilename+".pos");
                        if (inConfigFile.exists()) {
                            System.out.println("bringing configuration from v^2="+(v2+0.01*i));
                            ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
                            ConfigurationFileBinary.rescale(configFile, sim.box, 1/Math.sqrt(v2+0.01*i), Space3D.getInstance());
                            break;
                        }
                    }
                }
            }
        }
        
        if (!graphics && calcMu) {
            // we are calculating chemical potential via Widom insertion
            // mu ~= dA/dN ~= A(N+1) - A(N)
            // we introduce a finite size effect because the difference is centered about a system with slightly higher density
            // density => density * (N+1)/N
            // adjust nominal density so that
            // density * N/(N+0.5) => density * (N+1)/(N+0.5)
            // which is now centered about "density"
//            density = (density * numAtoms) / (numAtoms + 0.5);
            
            System.out.println("adjusting nominal density => "+density);
            throw new RuntimeException("fix me");
        }

        sim.potentialMasterList.getNeighborManager(sim.box).reset();

        if (!graphics) {
            long eqSteps = steps/10;
            if (eqSteps < 4000) {
                eqSteps = steps/4;
                if (eqSteps > 4000) eqSteps = 4000;
            }
            sim.ai.setMaxSteps(eqSteps);
            sim.getController().actionPerformed();
            sim.getController().reset();

            System.out.println("equilibration finished ("+eqSteps+" steps)");
        }
    	
        final MeterPotentialEnergy meterEnergyFast = new MeterPotentialEnergy(sim.potentialMasterList);
        meterEnergyFast.setBox(sim.box);
        long bs = steps/(fastInterval*nAccBlocks);
        final AccumulatorAverageFixed avgEnergyFast = new AccumulatorAverageFixed(bs == 0 ? 1 : bs);
        DataPumpListener energyPumpFast = new DataPumpListener(meterEnergyFast, avgEnergyFast, fastInterval);
        if (graphics) {
            sim.integrator.getEventManager().addListener(energyPumpFast);
        }

    	
    	MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMasterLong);
        meterEnergy.setBox(sim.box);
        final double uFac = meterEnergy.getDataAsScalar() - meterEnergyFast.getDataAsScalar();
        DataProcessor energyReweight = new DataProcessor() {
            protected final DataDoubleArray data = new DataDoubleArray(4);
            
            public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
                return null;
            }
            
            protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                dataInfo = new DataInfoDoubleArray("energy", Energy.DIMENSION, new int[]{4});
                return dataInfo;
            }

            protected IData processData(IData inputData) {
                double[] x = data.getData();
                x[0] = inputData.getValue(0);
                double uFast = meterEnergyFast.getDataAsScalar();

                double dx = x[0] - (uFast+uFac);
                double fac = 1;
                double w = Math.exp(-dx/temperature)/fac;
                x[1] = x[0] * w;
                x[2] = w;
                x[3] = uFast;
                return data;
            }
        };
        DataPumpListener energyPump = new DataPumpListener(meterEnergy, energyReweight, longInterval);
        if (graphics || !calcMu) {
            sim.integrator.getEventManager().addListener(energyPump);
        }
        bs = steps/(longInterval*nAccBlocks);
        final AccumulatorRatioAverageCovarianceFull avgEnergy = new AccumulatorRatioAverageCovarianceFull(bs == 0 ? 1 : bs);
        avgEnergy.setPushInterval(1);
        energyReweight.setDataSink(avgEnergy);

        // blocks should contain at least 25 samples (so they are mostly uncorrelated)
        long nLogSamples = steps/(longInterval*25);
        if (nLogSamples < 32) nLogSamples = 32;
        final int nDoubles = (int)(Math.log(nLogSamples)/Math.log(2));

        double rcMax = 0.494*L;
        double fac = 1.2;
        int nCutoffs = 1 + (int)(Math.log(rcMax/rcShort)/Math.log(fac)); 
        final double[] cutoffs = new double[nCutoffs];
        cutoffs[0] = rcShort;
        for (int i=1; i<cutoffs.length; i++) {
            cutoffs[i] = cutoffs[i-1]*1.2;
        }
        MeterPotentialEnergyCutoff meterEnergyCut = new MeterPotentialEnergyCutoff(sim.potentialMasterLongCut, sim.getSpace(), cutoffs);
        meterEnergyCut.setBox(sim.box);
        final double[] uFacCut = new double[cutoffs.length];
        IData uCut = meterEnergyCut.getData();
        double uFast0 = meterEnergyFast.getDataAsScalar();
        for (int i=0; i<cutoffs.length; i++) {
            uFacCut[i] = uCut.getValue(i) - uFast0;
        }
        DataProcessor energyReweightCut = new DataProcessor() {
            protected DataDoubleArray data;
            
            public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
                return null;
            }
            
            protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                data = new DataDoubleArray(inputDataInfo.getLength());
                dataInfo = new DataInfoDoubleArray("energy", Null.DIMENSION, new int[]{inputDataInfo.getLength()});
                return dataInfo;
            }
            
            protected IData processData(IData inputData) {
                double[] x = data.getData();
                double uFast = meterEnergyFast.getDataAsScalar();
                for (int i=0; i<x.length; i++) {
                    x[i] = Math.exp(-(inputData.getValue(i) - (uFast+uFacCut[i]))/temperature);
                }
                return data;
            }
        };
        DataPumpListener energyPumpCut = new DataPumpListener(meterEnergyCut, energyReweightCut, longInterval);
        if (graphics || !calcMu) {
            sim.integrator.getEventManager().addListener(energyPumpCut);
        }
        bs = steps/(longInterval*nAccBlocks);
        DataSplitter uCutSplitter = new DataSplitter();
        uCutSplitter.setDataSinkFactory(new IDataSinkFactory() {
            public IDataSink makeDataSink(int i) {
                AccumulatorAverageCollapsingLog acc = new AccumulatorAverageCollapsingLog(2);
                acc.setNumRawDataDoubles(nDoubles);
                return acc;
            }
        });
        energyReweightCut.setDataSink(uCutSplitter);
        

        final MeterPotentialEnergy meterEnergy2 = new MeterPotentialEnergy(sim.potentialMasterLong);
        meterEnergy2.setBox(sim.box);
        final ValueCache energyFastCache = new ValueCache(meterEnergyFast, sim.integrator);
        final ValueCache energyFullCache = new ValueCache(meterEnergy, sim.integrator);


        final MeterPressure meterPressureFast = new MeterPressure(sim.getSpace());
        meterPressureFast.setIntegrator(sim.integrator);
        bs = steps/(fastInterval*nAccBlocks);
        final AccumulatorAverageFixed avgPressureFast = new AccumulatorAverageFixed(bs == 0 ? 1 : bs);
        DataPumpListener pressurePumpFast = new DataPumpListener(meterPressureFast, avgPressureFast, fastInterval);
        if (graphics) {
            sim.integrator.getEventManager().addListener(pressurePumpFast);
        }

        
        MeterPressure meterPressure = new MeterPressure(sim.getSpace());
        meterPressure.setPotentialMaster(sim.potentialMasterLong);
        meterPressure.setTemperature(temperature);
        meterPressure.setBox(sim.box);
        DataProcessor pressureReweight = new DataProcessorReweight(temperature, energyFastCache, energyFullCache, uFac, Double.NaN, Double.NaN, sim.box, null);
        DataPumpListener pressurePump = new DataPumpListener(meterPressure, pressureReweight, longInterval);
        if (graphics || !calcMu) {
            sim.integrator.getEventManager().addListener(pressurePump);
        }
        bs = steps/(longInterval*nAccBlocks);
        final AccumulatorRatioAverageCovarianceFull avgPressure = new AccumulatorRatioAverageCovarianceFull(bs == 0 ? 1 : bs);
        avgPressure.setPushInterval(1);
        pressureReweight.setDataSink(avgPressure);

        final MeterPU meterPU = new MeterPU(sim.getSpace());
        meterPU.setBox(sim.box);
        meterPU.setIncludeLrc(false);
        meterPU.setPotentialMaster(sim.potentialMasterList);
        meterPU.setTemperature(temperature);
        
        bs = steps/(fastInterval*nAccBlocks);
        final AccumulatorRatioAverageCovarianceFull avgPU = new AccumulatorRatioAverageCovarianceFull(bs == 0 ? 1 : bs);
        DataPumpListener pressurePUFast = new DataPumpListener(meterPU, avgPU, fastInterval);
        if (graphics || !calcMu) {
            sim.integrator.getEventManager().addListener(pressurePUFast);
        }


        MeterPotentialEnergy meterEnergyWidomFast = new MeterPotentialEnergy(sim.potentialMasterList);
        MeterWidomInsertion meterWidomFast = new MeterWidomInsertion(sim.getSpace(), sim.getRandom());

        meterWidomFast.setEnergyMeter(meterEnergyWidomFast);
        meterWidomFast.setNInsert(numFastInsert);
        meterWidomFast.setSpecies(sim.species);
        meterWidomFast.setBox(sim.box);
        meterWidomFast.setTemperature(temperature);
        final AccumulatorAverageCollapsingLog avgLogWidomFast = new AccumulatorAverageCollapsingLog(2);
        // we are going to assume that each sample from widom is independent (allow full reblocking)
        bs = steps/(hybridInterval*nAccBlocks);
        final AccumulatorAverageFixed avgWidomFast = new AccumulatorAverageFixed(bs == 0 ? 1 : bs);
        avgWidomFast.setPushInterval(1);
        DataPumpListener pumpWidomFast = new DataPumpListener(meterWidomFast, graphics ? avgWidomFast : avgLogWidomFast, fastInterval);
        if (calcMu) {
            sim.integrator.getEventManager().addListener(pumpWidomFast);
        }

        bs = steps/(hybridInterval*10*nAccBlocks);
        
        MeterWidomInsertionCorrection meterWidomCorrection = new MeterWidomInsertionCorrection(sim.getSpace(), sim.getRandom());
        meterWidomCorrection.setEnergyMeter(meterEnergyWidomFast, meterEnergy2);
        meterWidomCorrection.setNInsert(numFullInsert);
        meterWidomCorrection.setSpecies(sim.species);
        meterWidomCorrection.setBox(sim.box);
        meterWidomCorrection.setTemperature(temperature);
        meterWidomCorrection.setEnergyFac(uFac);

        bs = steps/(longInterval*nAccBlocks);
        final AccumulatorRatioAverageCovarianceFull avgWidomCorrection = new AccumulatorRatioAverageCovarianceFull(bs == 0 ? 1 : bs);
        avgWidomCorrection.setPushInterval(1);
        DataPumpListener pumpWidomCorrection = new DataPumpListener(meterWidomCorrection, avgWidomCorrection, longInterval);
        if (calcMu) {
            sim.integrator.getEventManager().addListener(pumpWidomCorrection);
        }
        
        DataProcessor dpWidomCorrection = new DataProcessor() {
            protected final DataDouble data = new DataDouble();
            
            public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
                return null;
            }
            
            protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                dataInfo = new DataInfoDouble("mu", Null.DIMENSION);
                dataInfo.addTags(inputDataInfo.getTags());
                return dataInfo;
            }
            
            protected IData processData(IData inputData) {
                double avgFast = avgWidomFast.getData(avgWidomFast.AVERAGE).getValue(0);
                data.x = avgFast;
                if (avgFast > 0) {
                    double avgFull = inputData.getValue(0);
                    double avgFastFull = inputData.getValue(1);
                    double avgFullCorrection = inputData.getValue(2);
//                    System.out.println(String.format("%15.5e %15.5e %15.5e %15.5e", avgFast, avgFull, avgFastFull, avgFullCorrection));
                    data.x *= avgFull/(avgFastFull*avgFullCorrection);
//                    System.out.println(avgFull/(avgFastFull*avgFullCorrection)+" "+inputData.getValue(3)/inputData.getValue(1));
                    if (data.x == Double.POSITIVE_INFINITY) {
                        // just preliminary data, or we need to adjust energyFac
                        data.x = 0;
                    }
                }
                return data;
            }
        };
        
        avgWidomCorrection.addDataSink(dpWidomCorrection, new StatType[]{avgWidomCorrection.AVERAGE});
        Function muFunction = new Function() {
            public double f(double x) {
                if (x==0) return Double.NaN;
                return -temperature*Math.log(x);
            }
        };
        DataProcessorFunction muFast = new DataProcessorFunction(muFunction);
        avgWidomFast.addDataSink(muFast, new StatType[]{avgWidomFast.AVERAGE});
        DataProcessorFunction muCorrected = new DataProcessorFunction(muFunction);
        dpWidomCorrection.setDataSink(muCorrected);

        MeterWidomInsertionP meterWidomDFast = new MeterWidomInsertionP(sim.getSpace(), sim.getRandom());
        meterWidomDFast.setEnergyMeter(meterEnergyWidomFast);
        meterWidomDFast.setNInsert(numFastInsert);
        meterWidomDFast.setSpecies(sim.species);
        meterWidomDFast.setBox(sim.box);
        meterWidomDFast.setTemperature(temperature);
        meterWidomDFast.setPressure(Double.NaN);
        bs = steps/(fastInterval*nAccBlocks);
        final AccumulatorAverageFixed avgWidomDFast = new AccumulatorAverageFixed(bs == 0 ? 1 : bs);
        avgWidomDFast.setPushInterval(1);
        DataPumpListener pumpWidomDFast = new DataPumpListener(meterWidomDFast, avgWidomDFast, fastInterval);
        if (calcMu) {
            sim.integrator.getEventManager().addListener(pumpWidomDFast);
        }

        Function muDFunction = new Function() {
            public double f(double x) {
                if (x==0) return Double.NaN;
                return temperature*Math.log(x);
            }
        };
        DataProcessorFunction muDFast = new DataProcessorFunction(muFunction);
        avgWidomDFast.addDataSink(muDFast, new StatType[]{avgWidomFast.AVERAGE});

        
    	if (graphics) {
            final String APP_NAME = "LjMd3D";
        	final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());
    
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.getController().getDataStreamPumps().add(energyPump);
            simGraphic.getController().getDataStreamPumps().add(energyPumpFast);
            simGraphic.getController().getDataStreamPumps().add(pressurePump);
            simGraphic.getController().getDataStreamPumps().add(pressurePumpFast);
            simGraphic.getController().getDataStreamPumps().add(pumpWidomCorrection);
            simGraphic.getController().getDataStreamPumps().add(pumpWidomFast);
            simGraphic.getController().getDataStreamPumps().add(pumpWidomDFast);

            simGraphic.makeAndDisplayFrame(APP_NAME);
    
//            DisplayTextBoxesCAE displayEnergy = new DisplayTextBoxesCAE();
//            displayEnergy.setAccumulator(avgEnergy);
//            simGraphic.add(displayEnergy);
//            DisplayTextBoxesCAE displayEnergyFast = new DisplayTextBoxesCAE();
//            displayEnergyFast.setAccumulator(avgEnergyFast);
//            simGraphic.add(displayEnergyFast);
            
            AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSourceCountTime timeDataSource = new DataSourceCountTime(sim.integrator);
            energyHistory.setTimeDataSource(timeDataSource);
            energyHistory.setPushInterval(1);
            DataFork energyFork = new DataFork();
            energyPump.setDataSink(energyFork);
            energyFork.addDataSink(energyReweight);
            energyFork.addDataSink(energyHistory);
            
            AccumulatorHistory energyCorrectedAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            energyCorrectedAvgHistory.setTimeDataSource(timeDataSource);
            energyCorrectedAvgHistory.setPushInterval(1);
            DataProcessor dpAvgEnergyCorrected = new DataProcessor() {
                protected final DataDouble data = new DataDouble();
                
                public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {return null;}
                protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                    dataInfo = new DataInfoDouble("energy", Energy.DIMENSION);
                    dataInfo.addTags(inputDataInfo.getTags());
                    return dataInfo;
                }
                protected IData processData(IData inputData) {
                    IData avg = ((DataGroup)inputData).getData(avgEnergy.AVERAGE.index);
                    data.x = avg.getValue(1)/avg.getValue(2);
                    return data;
                }
            };
            avgEnergy.addDataSink(dpAvgEnergyCorrected);
            dpAvgEnergyCorrected.setDataSink(energyCorrectedAvgHistory);

            AccumulatorHistory energyCorrected2AvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            energyCorrected2AvgHistory.setTimeDataSource(timeDataSource);
            energyCorrected2AvgHistory.setPushInterval(1);
            DataProcessor dpAvgEnergyCorrected2 = new DataProcessor() {
                protected final DataDouble data = new DataDouble();
                
                public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {return null;}
                protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                    dataInfo = new DataInfoDouble("energy", Energy.DIMENSION);
                    dataInfo.addTags(inputDataInfo.getTags());
                    return dataInfo;
                }
                protected IData processData(IData inputData) {
                    IData avg = ((DataGroup)inputData).getData(avgEnergy.AVERAGE.index);
                    data.x = avg.getValue(1)/(avg.getValue(3)*avg.getValue(2)) * avgEnergyFast.getData(avgEnergyFast.AVERAGE).getValue(0);
                    return data;
                }
            };
            avgEnergy.addDataSink(dpAvgEnergyCorrected2);
            dpAvgEnergyCorrected2.setDataSink(energyCorrected2AvgHistory);

            AccumulatorHistory energyAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            energyAvgHistory.setTimeDataSource(timeDataSource);
            energyAvgHistory.setPushInterval(1);
            DataProcessor dpAvgEnergy = new DataProcessor() {
                protected final DataDouble data = new DataDouble();
                
                public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {return null;}
                protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
                    dataInfo = new DataInfoDouble("energy", Energy.DIMENSION);
                    dataInfo.addTags(inputDataInfo.getTags());
                    return dataInfo;
                }
                protected IData processData(IData inputData) {
                    IData avg = ((DataGroup)inputData).getData(avgEnergy.AVERAGE.index);
                    data.x = avg.getValue(0);
                    return data;
                }
            };
            avgEnergy.addDataSink(dpAvgEnergy);
            dpAvgEnergy.setDataSink(energyAvgHistory);
            
            AccumulatorHistory energyHistoryFast = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyHistoryFast.setTimeDataSource(timeDataSource);
            energyHistoryFast.setPushInterval(1);
            DataFork energyForkFast = new DataFork();
            energyPumpFast.setDataSink(energyForkFast);
            energyForkFast.addDataSink(avgEnergyFast);
            energyForkFast.addDataSink(energyHistoryFast);

            DisplayPlot energyPlot = new DisplayPlot();
            
            energyHistoryFast.setDataSink(energyPlot.getDataSet().makeDataSink());
            energyHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
            energyAvgHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
            energyCorrectedAvgHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
            energyCorrected2AvgHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
            energyPlot.setLabel("energy");
            energyPlot.setLegend(new DataTag[]{energyHistoryFast.getTag()}, "fast energy");
            energyPlot.setLegend(new DataTag[]{energyHistory.getTag()}, "energy");
            energyPlot.setLegend(new DataTag[]{energyAvgHistory.getTag()}, "avg energy");
            energyPlot.setLegend(new DataTag[]{energyCorrectedAvgHistory.getTag()}, "avg corrected energy");
            energyPlot.setLegend(new DataTag[]{energyCorrected2AvgHistory.getTag()}, "avg corrected2 energy");
            energyPlot.setUnit(new SimpleUnit(Energy.DIMENSION, numAtoms, "energy", "", false));
            simGraphic.add(energyPlot);
            
            
            AccumulatorHistory pressureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            pressureHistory.setTimeDataSource(timeDataSource);
            pressureHistory.setPushInterval(1);
            DataFork pressureFork = new DataFork();
            pressurePump.setDataSink(pressureFork);
            pressureFork.addDataSink(pressureReweight);
            pressureFork.addDataSink(pressureHistory);
            
            AccumulatorHistory pressureCorrectedAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            pressureCorrectedAvgHistory.setTimeDataSource(timeDataSource);
            pressureCorrectedAvgHistory.setPushInterval(1);
            DataProcessor dpAvgPressureCorrected = new DataProcessorCorrection(1);
            avgPressure.addDataSink(dpAvgPressureCorrected);
            dpAvgPressureCorrected.setDataSink(pressureCorrectedAvgHistory);

            AccumulatorHistory pressureAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            pressureAvgHistory.setTimeDataSource(timeDataSource);
            pressureAvgHistory.setPushInterval(1);
            AccumulatorAverageFixed avgRawPressure = new AccumulatorAverageFixed(bs);
            pressureFork.addDataSink(avgRawPressure);
            avgRawPressure.setPushInterval(1);
            avgRawPressure.addDataSink(pressureAvgHistory, new StatType[]{avgRawPressure.AVERAGE});
            
            AccumulatorHistory pressureHistoryFast = new AccumulatorHistory(new HistoryCollapsingAverage());
            pressureHistoryFast.setTimeDataSource(timeDataSource);
            pressureHistoryFast.setPushInterval(1);
            DataFork pressureForkFast = new DataFork();
            pressurePumpFast.setDataSink(pressureForkFast);
            pressureForkFast.addDataSink(avgPressureFast);
            pressureForkFast.addDataSink(pressureHistoryFast);

            DisplayPlot pressurePlot = new DisplayPlot();
            
            pressureHistoryFast.setDataSink(pressurePlot.getDataSet().makeDataSink());
            pressureHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
            pressureAvgHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
            pressureCorrectedAvgHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
            pressurePlot.setLabel("pressure");
            pressurePlot.setLegend(new DataTag[]{pressureHistoryFast.getTag()}, "fast pressure");
            pressurePlot.setLegend(new DataTag[]{pressureHistory.getTag()}, "pressure");
            pressurePlot.setLegend(new DataTag[]{pressureAvgHistory.getTag()}, "avg pressure");
            pressurePlot.setLegend(new DataTag[]{pressureCorrectedAvgHistory.getTag()}, "avg corrected pressure");
            simGraphic.add(pressurePlot);


            if (calcMu) {
                AccumulatorHistory widomHistoryFast = new AccumulatorHistory(new HistoryCollapsingDiscard());
                widomHistoryFast.setTimeDataSource(timeDataSource);
                widomHistoryFast.setPushInterval(1);
                muFast.setDataSink(widomHistoryFast);
    
                AccumulatorHistory widomCorrectionHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
                widomCorrectionHistory.setTimeDataSource(timeDataSource);
                widomCorrectionHistory.setPushInterval(1);
                muCorrected.setDataSink(widomCorrectionHistory);
    
                AccumulatorHistory widomDHistoryFast = new AccumulatorHistory(new HistoryCollapsingDiscard());
                widomDHistoryFast.setTimeDataSource(timeDataSource);
                widomDHistoryFast.setPushInterval(1);
                muDFast.setDataSink(widomDHistoryFast);
    
    
                DisplayPlot widomPlot = new DisplayPlot();
                widomHistoryFast.setDataSink(widomPlot.getDataSet().makeDataSink());
                widomCorrectionHistory.setDataSink(widomPlot.getDataSet().makeDataSink());
                widomDHistoryFast.setDataSink(widomPlot.getDataSet().makeDataSink());
                widomPlot.setLabel("widom");
                widomPlot.setLegend(new DataTag[]{widomHistoryFast.getTag()}, "widom fast");
                widomPlot.setLegend(new DataTag[]{widomCorrectionHistory.getTag()}, "widom corrected");
                widomPlot.setLegend(new DataTag[]{widomDHistoryFast.getTag()}, "widom fast deletion");
                simGraphic.add(widomPlot);
            }
            
            
            return;
    	}

    	long t1 = System.currentTimeMillis();
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        System.out.println();
        
        WriteConfigurationBinary writeConfig = new WriteConfigurationBinary(sim.getSpace());
        writeConfig.setFileName(configFilename+".pos");
        writeConfig.setBox(sim.box);
        writeConfig.actionPerformed();
        
        System.out.println("hybrid acceptance: "+sim.integrator.getHybridAcceptance());
        System.out.println();


        if (!calcMu) {
            double uFastAvg = avgEnergyFast.getData(avgEnergyFast.AVERAGE).getValue(0);
            double uFastErr = avgEnergyFast.getData(avgEnergyFast.ERROR).getValue(0);
            double uFastCor = avgEnergyFast.getData(avgEnergyFast.BLOCK_CORRELATION).getValue(0);
    
//            System.out.println(String.format("potential energy (fast): %20.15e   %10.4e   %4.2f", uFastAvg/numAtoms, uFastErr/numAtoms, uFastCor));

            IData uAvgData = avgEnergy.getData(avgEnergy.AVERAGE);
            IData uRatioData = avgEnergy.getData(avgEnergy.RATIO);
            IData uErrData = avgEnergy.getData(avgEnergy.ERROR);
            IData uRatioErrData = avgEnergy.getData(avgEnergy.RATIO_ERROR);
            IData uCorData = avgEnergy.getData(avgEnergy.BLOCK_CORRELATION);
            IData uCovData = avgEnergy.getData(avgEnergy.BLOCK_COVARIANCE);
    
            double uFullAvg = uAvgData.getValue(0);
            double uFullErr = uErrData.getValue(0);
            double uFullCor = uCorData.getValue(0);
    
            System.out.println(avgEnergy.getBlockCount()+" energy blocks");
            System.out.println(String.format("potential energy (full): %20.15e   %10.4e   %4.2f", uFullAvg/numAtoms, uFullErr/numAtoms, uFullCor));

            double uReweightedAvg = uRatioData.getValue(4*1+2);
            double uReweightedErr = uRatioErrData.getValue(4*1+2);
            double uReweightedCor = uCorData.getValue(1);
            if (uCorData.getValue(2) > uReweightedCor) uReweightedCor = uCorData.getValue(2);
    
            System.out.println(String.format("potential energy (reweighted): %20.15e   %10.4e   %4.2f", uReweightedAvg/numAtoms, uReweightedErr/numAtoms, uReweightedCor));

        
            if (longInterval > hybridInterval) {
                double uReweighted2Avg = uAvgData.getValue(1)/(uAvgData.getValue(2)*uAvgData.getValue(3)) * uFastAvg;
                double uReweighted2Err = Math.pow(uErrData.getValue(1)/uAvgData.getValue(1), 2) + Math.pow(uErrData.getValue(2)/uAvgData.getValue(2), 2) +
                                                                    Math.pow(uErrData.getValue(3)/uAvgData.getValue(3), 2) + Math.pow(uFastErr/uFastAvg, 2);
                long blockCount = avgEnergy.getBlockCount();
                uReweighted2Err += -2*uCovData.getValue(4*1+2)/((blockCount-1)*uAvgData.getValue(1)*uAvgData.getValue(2));
                uReweighted2Err += -2*uCovData.getValue(4*1+3)/((blockCount-1)*uAvgData.getValue(1)*uAvgData.getValue(3));
                uReweighted2Err +=  2*uCovData.getValue(4*2+3)/((blockCount-1)*uAvgData.getValue(2)*uAvgData.getValue(3));
                uReweighted2Err = Math.abs(uReweighted2Avg) * Math.sqrt(uReweighted2Err);
                        
                double uReweighted2Cor = uReweightedCor;
                if (uCorData.getValue(3) > uReweighted2Cor) uReweighted2Cor = uCorData.getValue(3);
                if (uFastCor > uReweighted2Cor) uReweighted2Cor = uFastCor;
        
                System.out.println(String.format("potential energy (reweighted fast): %20.15e   %10.4e   %4.2f", uReweighted2Avg/numAtoms, uReweighted2Err/numAtoms, uReweighted2Cor));
            }

            double pFastAvg = avgPressureFast.getData(avgPressureFast.AVERAGE).getValue(0);
            double pFastErr = avgPressureFast.getData(avgPressureFast.ERROR).getValue(0);
            double pFastCor = avgPressureFast.getData(avgPressureFast.BLOCK_CORRELATION).getValue(0);

//            System.out.println(String.format("pressure (fast): %20.15e   %10.4e   %4.2f", pFastAvg, pFastErr, pFastCor));

            IData pAvgData = avgPressure.getData(avgPressure.AVERAGE);
            IData pRatioData = avgPressure.getData(avgPressure.RATIO);
            IData pErrData = avgPressure.getData(avgPressure.ERROR);
            IData pRatioErrData = avgPressure.getData(avgPressure.RATIO_ERROR);
            IData pCorData = avgPressure.getData(avgPressure.BLOCK_CORRELATION);
            IData pCovData = avgPressure.getData(avgPressure.BLOCK_COVARIANCE);

            double pFullAvg = pAvgData.getValue(0);
            double pFullErr = pErrData.getValue(0);
            double pFullCor = pCorData.getValue(0);

            System.out.println("\n"+avgPressure.getBlockCount()+" pressure blocks");
            System.out.println(String.format("pressure (full): %20.15e   %10.4e   %4.2f", pFullAvg, pFullErr, pFullCor));

            double pReweightedAvg = pRatioData.getValue(1);
            double pReweightedErr = pRatioErrData.getValue(1);
            double pReweightedCor = pCorData.getValue(1);
//                if (pCorData.getValue(2) > pReweightedCor) pReweightedCor = pCorData.getValue(2);

            System.out.println(String.format("pressure (reweighted): %20.15e   %10.4e   %4.2f", pReweightedAvg, pReweightedErr, pReweightedCor));

            if (longInterval > hybridInterval) {
                double pReweighted2Avg = pAvgData.getValue(1)/(pAvgData.getValue(2)*pAvgData.getValue(3)) * pFastAvg;
                double pReweighted2Err = Math.pow(pErrData.getValue(1)/pAvgData.getValue(1), 2) + Math.pow(pErrData.getValue(2)/pAvgData.getValue(2), 2) +
                                                                    Math.pow(pErrData.getValue(3)/pAvgData.getValue(3), 2) + Math.pow(pFastErr/pFastAvg, 2);
                long blockCount = avgEnergy.getBlockCount();
                pReweighted2Err += -2*pCovData.getValue(4*1+2)/((blockCount-1)*pAvgData.getValue(1)*pAvgData.getValue(2));
                pReweighted2Err += -2*pCovData.getValue(4*1+3)/((blockCount-1)*pAvgData.getValue(1)*pAvgData.getValue(3));
                pReweighted2Err +=  2*pCovData.getValue(4*2+3)/((blockCount-1)*pAvgData.getValue(2)*pAvgData.getValue(3));
                pReweighted2Err = Math.abs(pReweighted2Avg) * Math.sqrt(pReweighted2Err);
                        
                double pReweighted2Cor = pReweightedCor;
                if (pCorData.getValue(3) > pReweighted2Cor) pReweighted2Cor = pCorData.getValue(3);
                if (pFastCor > pReweighted2Cor) pReweighted2Cor = pFastCor;

                System.out.println(String.format("pressure (reweighted fast): %20.15e   %10.4e   %4.2f", pReweighted2Avg, pReweighted2Err, pReweighted2Cor));
            }

            System.out.println();
            
            for (int i=0; i<cutoffs.length; i++) {
                AccumulatorAverageCollapsingLog accFECut = (AccumulatorAverageCollapsingLog)uCutSplitter.getDataSink(i);
                IData feAvgCutData = accFECut.getAverages();
                int feCutLength = feAvgCutData.getLength();
                double feAvgCut = -temperature*Math.log(feAvgCutData.getValue(feCutLength-1));
                IData feLogErrCutData = accFECut.getStdevLog();
                double feErrCut = temperature*feLogErrCutData.getValue(feCutLength-1);
                double feCorCut = accFECut.getRawLogDataCorrelation();
                
                P2SoftSphericalTruncated p2t = new P2SoftSphericalTruncated(sim.getSpace(), sim.potential, cutoffs[i]);
                p2t.setBox(sim.box);
                Potential0Lrc plrc = p2t.makeLrcPotential(new IAtomType[]{sim.species.getAtomType(0), sim.species.getAtomType(0)});
                plrc.setBox(sim.box);
                double ulrc = plrc.energy(null);

                System.out.println(String.format("rc: %d  A-Afast: %20.15e   %10.4e  %4.2f (%d samples)", i, (ulrc + uFacCut[i] + feAvgCut)/numAtoms, feErrCut/numAtoms, feCorCut, accFECut.getNumRawData()));
            }
            System.out.println();

            IData puAvg = avgPU.getData(avgPU.AVERAGE);
            IData puErr = avgPU.getData(avgPU.ERROR);
            IData puCor = avgPU.getData(avgPU.BLOCK_CORRELATION);
            IData puCov = avgPU.getData(avgPU.COVARIANCE);
            System.out.println(String.format("U:     % 20.15e   %10.4e   %4.2f", puAvg.getValue(0), puErr.getValue(0), puCor.getValue(0)));
            double corPU = puCov.getValue(1)/Math.sqrt(puCov.getValue(0)*puCov.getValue(5));
            System.out.println(String.format("P:     % 20.15e   %10.4e   %4.2f  %6.4f", puAvg.getValue(1), puErr.getValue(1), puCor.getValue(1), corPU));
            System.out.println();
            System.out.println(String.format("DADy:  % 20.15e   %10.4e   %4.2f", puAvg.getValue(2), puErr.getValue(2), puCor.getValue(2)));
            double cor = puCov.getValue(11)/Math.sqrt(puCov.getValue(10)*puCov.getValue(15));
            System.out.println(String.format("DADv2: % 20.15e   %10.4e   %4.2f  %6.4f", puAvg.getValue(3), puErr.getValue(3), puCor.getValue(3), -cor));

        }
        else {
//            System.out.println(avgLogWidomFast.getCount()+" fast widom samples");
    
            IData muFastAvgData = avgLogWidomFast.getAverages();
            int muLength = muFastAvgData.getLength();
            double muFastAvg = -temperature*Math.log(muFastAvgData.getValue(muLength-1));
            IData muFastLogErrData = avgLogWidomFast.getStdevLog();
            double muFastErr = temperature*muFastLogErrData.getValue(muLength-1);
            double muCor = avgLogWidomFast.getRawLogDataCorrelation();
    
            System.out.println(String.format("chemical potential (fast): %20.15e   %10.4e  %4.2f  (%d samples)", muFastAvg, muFastErr, muCor, avgLogWidomFast.getNumRawData()));
        }
/*
        if (false) {
        IData muAvgData = avgWidomCorrection.getData(avgWidomCorrection.AVERAGE);
        IData muErrData = avgWidomCorrection.getData(avgWidomCorrection.ERROR);
        IData muCorData = avgWidomCorrection.getData(avgWidomCorrection.BLOCK_CORRELATION);
        IData muRatioData = avgWidomCorrection.getData(avgWidomCorrection.RATIO);
        IData muRatioErrData = avgWidomCorrection.getData(avgWidomCorrection.RATIO_ERROR);
        IData muCovData = avgWidomCorrection.getData(avgWidomCorrection.BLOCK_COVARIANCE);
        
        double muFullAvg = muAvgData.getValue(3);
        double muFullErr = muErrData.getValue(3);
        double muFullCor = muCorData.getValue(3);

        System.out.println(avgWidomCorrection.getBlockCount()+" full widom blocks");

        System.out.println(String.format("chemical potential (full): %20.15e   %10.4e   %4.2f", -Math.log(muFullAvg), muFullErr/muFullAvg, muFullCor));

        double muReweightedAvg = muRatioData.getValue(4*0+2);
        double muReweightedErr = muRatioErrData.getValue(4*0+2);
        double muReweightedCor = muCorData.getValue(0);
        if (muReweightedCor < muCorData.getValue(2)) muReweightedCor = muCorData.getValue(2);

        System.out.println(String.format("chemical potential (reweighted): %20.15e   %10.4e   %4.2f", -Math.log(muReweightedAvg), muReweightedErr/muReweightedAvg, muReweightedCor));

        double muReweighted2Avg = muReweightedAvg / muAvgData.getValue(1) * muFastAvg;
        double muReweighted2Err = Math.pow(muErrData.getValue(0)/muAvgData.getValue(0), 2) + Math.pow(muErrData.getValue(1)/muAvgData.getValue(1), 2) 
                                + Math.pow(muErrData.getValue(2)/muAvgData.getValue(2), 2);
        long blockCount = avgWidomCorrection.getBlockCount();
        muReweighted2Err += -2*muCovData.getValue(4*0+1)/((blockCount-1)*muAvgData.getValue(0)*muAvgData.getValue(1));
        muReweighted2Err += -2*muCovData.getValue(4*0+2)/((blockCount-1)*muAvgData.getValue(0)*muAvgData.getValue(2));
        muReweighted2Err +=  2*muCovData.getValue(4*1+2)/((blockCount-1)*muAvgData.getValue(1)*muAvgData.getValue(2));
        muReweighted2Err += Math.pow(muFastErr/muFastAvg, 2);
        muReweighted2Err = muReweighted2Avg * Math.sqrt(muReweighted2Err);
        double muReweighted2Cor = muReweightedCor;
        if (muReweightedCor < muCorData.getValue(1)) muReweightedCor = muCorData.getValue(1);

        System.out.println(String.format("chemical potential (reweighted fast): %20.15e   %10.4e   %4.2f", -Math.log(muReweighted2Avg), muReweighted2Err/muReweighted2Avg, muReweighted2Cor));
        }
        */

        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }
    
    public static class DataProcessorCorrection extends DataProcessor {
        protected DataDoubleArray data;
        protected final int nMu;
        
        public DataProcessorCorrection(int nMu) {
            this.nMu = nMu;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {return null;}

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataInfoDoubleArray("foo", Null.DIMENSION, new int[]{((DataInfoGroup)inputDataInfo).getSubDataInfo(AccumulatorAverageCovariance.AVERAGE.index).getLength()-nMu});
            data = new DataDoubleArray(dataInfo.getLength());
            return dataInfo;
        }

        protected IData processData(IData inputData) {
            IData avg = ((DataGroup)inputData).getData(AccumulatorAverageCovariance.AVERAGE.index);
            double[] x = data.getData();
            int nValues = avg.getLength()/nMu;
            for (int i=0; i<nMu; i++) {
                double wAvg = avg.getValue((i+1)*nValues-1);
                for (int j=0; j<nValues-1; j++) {
                    x[j+i*(nValues-1)] = avg.getValue(j+i*nValues)/wAvg;
                }
            }
            return data;
        }
    }

    public static class ValueCache {
        protected long lastStep = -1;
        protected double lastValue;
        protected final DataSourceScalar dss;
        protected final IIntegrator integrator;
        public ValueCache(DataSourceScalar dss, IIntegrator integrator) {
            this.dss = dss;
            this.integrator = integrator;
        }
        public double getValue() {
            if (integrator.getStepCount() != lastStep) {
                lastStep = integrator.getStepCount();
                lastValue = dss.getDataAsScalar();
            }
            return lastValue;
        }
    }
    
    public static class DataProcessorReweight extends DataProcessor {
        private final double temperature;
        private final ValueCache energyFastCache;
        private final ValueCache energyFullCache;
        private final double uFac;
        protected DataDoubleArray data;
        protected final double deltaMu, deltaP;
        protected final IBox box;
        protected final int numMolecules0;
        protected final double v0;
        protected final IFunction vBias;

        public DataProcessorReweight(double temperature,
                ValueCache energyFastCache, ValueCache energyFullCache,
                double uFac, double deltaMu, double deltaP, IBox box, IFunction vBias) {
            this.temperature = temperature;
            this.energyFastCache = energyFastCache;
            this.energyFullCache = energyFullCache;
            this.uFac = uFac;
            this.deltaMu = deltaMu;
            this.deltaP = deltaP;
            this.box = box;
            this.v0 = box.getBoundary().volume();
            this.numMolecules0 = box.getMoleculeList().getMoleculeCount();
            this.vBias = vBias;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
            return null;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataInfoDoubleArray("muX", Null.DIMENSION, new int[]{(inputDataInfo.getLength()+1)});
            data = new DataDoubleArray(dataInfo.getLength());
            return dataInfo;
        }

        protected IData processData(IData inputData) {
            double uFast = energyFastCache.getValue();
            double uFull = energyFullCache.getValue();
            int numMolecules = box.getMoleculeList().getMoleculeCount();
            double[] x = data.getData();
            double dx = uFull - (uFast+uFac);
            double fac = 1;
            if (!Double.isNaN(deltaMu)) {
                dx -= deltaMu*(numMolecules-numMolecules0);
            }
            else if (!Double.isNaN(deltaP)) {
                dx += deltaP*(box.getBoundary().volume() - v0);
                fac = vBias.f(box.getBoundary().volume());
            }
            double w = Math.exp(-dx/temperature)/fac;
            for (int i=0; i<inputData.getLength(); i++) {
                x[i] = inputData.getValue(i)*w;
            }
            x[inputData.getLength()] = w;
            if (data.isNaN()) {
                throw new RuntimeException("oops");
            }
            return data;
        }
    }

    public static class LjMd3DParams extends ParameterBase {
        public int numAtoms = 2048;
        public double y = 1.3;
        public double v2 = 0.8;
        public double tStepFac = 1;
        public long steps = 40000;
        public int hybridInterval = 20;
        public double rcShort = 2.5;
        public boolean graphics = false;
        public int nAccBlocks = 100;
        public boolean calcMu = false;
    }
}
