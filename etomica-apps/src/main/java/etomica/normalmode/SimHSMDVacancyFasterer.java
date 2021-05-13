/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHardFasterer;
import etomica.graphics.*;
import etomica.integrator.*;
import etomica.integrator.mcmove.MCMoveIDBiasActionFasterer;
import etomica.integrator.mcmove.MCMoveInsertDeleteLatticeVacancyFasterer;
import etomica.integrator.mcmove.MCMoveOverlapListenerFasterer;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.modifier.Modifier;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.normalmode.DataSourceMuRootFasterer.DataSourceMuRootVacancyConcentration;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.statmech.HardSphereSolid;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimHSMDVacancyFasterer extends Simulation {

    public final PotentialComputePair potentialMasterList;
    public final NeighborListManagerFastererHard neighborManager;
    public IntegratorHardMDMCFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2HardGeneric potential;
    public IntegratorMCFasterer integratorMC;
    public MCMoveVolume mcMoveVolume;
    public MCMoveInsertDeleteLatticeVacancyFasterer mcMoveID;


    public SimHSMDVacancyFasterer(final int numAtoms, double density, double tStep, int hybridInterval, final int numV, final double mu) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        PrimitiveCubic primitive = new PrimitiveCubic(space, n * L);
        int[] nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        BasisBigCell basis = new BasisBigCell(space, basisFCC, nCells);


        double nbrRange = 1.7;
        neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, nbrRange, BondingInfo.noBonding());
        potentialMasterList = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        double sigma = 1.0;
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        if (nbrRange > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("rcShort is too large");
        }

        double nbr1 = L / Math.sqrt(2);
        double y = 1.25 * nbr1; //nbr1+(L-nbr1)*0.6+0.06;

        potential = P2HardSphere.makePotential(y);
        AtomType leafType = species.getLeafType();

        potentialMasterList.setPairPotential(leafType, leafType, potential);
        integrator = new IntegratorHardMDMCFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMasterList), neighborManager, random, tStep, 1.0, box);
        integrator.setIsothermal(true);
        integrator.setThermostatNoDrift(true);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.HYBRID_MC);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        integrator.setThermostatInterval(hybridInterval);

        // we just needed this to initialize coordinates, not keep track of lattice positions of atoms
        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        NeighborCellManagerFasterer neighborCellManager = new NeighborCellManagerFasterer(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        PotentialComputePair potentialComputeMC = new PotentialComputePair(getSpeciesManager(), box, neighborCellManager, potentialMasterList.getPairPotentials());
        neighborCellManager.setPotentialRange(nbrRange);
        integratorMC = new IntegratorMCFasterer(potentialComputeMC, random, 1, box);
        mcMoveID = new MCMoveInsertDeleteLatticeVacancyFasterer(potentialComputeMC, neighborCellManager, random, space, integrator, y, numAtoms, numV);
        double x = (nbr1 - 1) / 4.0;
        mcMoveID.setMaxInsertDistance(x);
        mcMoveID.makeFccVectors(nbr1);
        mcMoveID.setMu(mu);
        mcMoveID.setSpecies(species);
        integratorMC.getMoveManager().addMCMove(mcMoveID);
        integratorMC.getMoveEventManager().addListener(mcMoveID);
        integrator.setIntegratorMC(integratorMC, numAtoms);
        potentialMasterList.init();

        // we set the potential range high so that the MC move could find its neighbors.
        // now set potential range back to its appropriate value.
        potential.setCollisionDiameter(0, sigma);
    }

    public static void main(String[] args) {

        HSMDVParams params = new HSMDVParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.numAtoms = 4000;
            params.steps = 1000000;
            params.density = 1.0;
            params.numV = 10;
        }

        final int numAtoms = params.numAtoms;
        final double density = params.density;
        long steps = params.steps;
        double tStep = params.tStep;
        double beta = (1 - density/Math.sqrt(2));
        tStep *= beta/numAtoms;
        boolean graphics = params.graphics;
        boolean fluid = params.fluid;
        int hi = params.hybridInterval;
        if (hi == -1) {
            // this keeps computational cost of insertions balanced with MD
            hi = numAtoms/10;
        }
        int hybridInterval = hi;
        int numV = params.numV;
        double mu = params.mu;
        double Alat = params.Alat;
        System.out.println("Running N="+params.numAtoms+" at rho="+density);
        if (Double.isNaN(mu)) {
            if (!fluid) {
                // this should get us pretty close
                Alat = HardSphereSolid.Asolid(density, 100000);
                System.out.println("lattice pressure "+HardSphereSolid.zSolid(density)*density);
                mu = (Alat + HardSphereSolid.zSolid(density));
                System.out.println("mu "+mu);
            }
            else {
                Alat = 5;
                mu = 8;
            }
        }
        else if (Double.isNaN(Alat)) {
            if (!fluid) {
                // this should get us pretty close
                Alat = HardSphereSolid.Asolid(density, 100000);
            }
            else {
                Alat = 5;
            }
        }
        double daDef = params.daDef;
        if (Double.isNaN(daDef)) {
            // use a correlation.  this should be with ~1
            double x = (1.0-density/Math.sqrt(2));
            daDef = 32.4 - 37*Math.pow(x, -0.065);
            System.out.println("estimated defect free energy: "+daDef);
        }
        else {
            System.out.println("using defect free energy: "+daDef);
        }
   
        if (!graphics) {
    	    System.out.println("Running HS MD with N="+numAtoms+" at density="+density);
    	    System.out.println("time step: "+tStep);
    	    System.out.println(steps+" steps ("+(steps*tStep)+" time units)");
    	    System.out.println("hybrid MC interval: "+hybridInterval);
    	}

        int fac = 1;
        // tStep needs to be less than 0.01 for happy neighbor updating 
        double tStepMin = fluid ? 0.002 : 0.005;
        while (tStep < tStepMin/10 && steps >= 1000*10) {
            fac *= 10;
            tStep *= 10;
            steps /= 10;
        }
        while (tStep < tStepMin/2 && steps > 1000*2 && ((steps/2000)*2000 == steps)) {
            fac *= 2;
            tStep *= 2;
            steps /= 2;
        }
        hybridInterval /= fac;
        if (fac != 1) {
            System.out.println(" => "+steps+" x "+tStep);
        }

        final SimHSMDVacancyFasterer sim = new SimHSMDVacancyFasterer(numAtoms, density, tStep, hybridInterval, numV, mu);


        final int biasInterval = 1000;

        final MCMoveOverlapListenerFasterer mcMoveOverlapMeter = new MCMoveOverlapListenerFasterer(sim.mcMoveID, 11, daDef, numAtoms, 2);
        mcMoveOverlapMeter.setTemperature(1);
        sim.integratorMC.getMoveEventManager().addListener(mcMoveOverlapMeter);

        final MeterPressureHardFasterer meterP = new MeterPressureHardFasterer(sim.integrator);
        DataDistributer.Indexer indexer = new DataDistributer.Indexer() {
            public int getIndex() {
                return numAtoms-sim.box.getNMolecules(sim.species);
            }
        };
        final DataDistributer pSplitter = new DataDistributer(indexer, new IDataSinkFactory() {
            public IDataSink makeDataSink(int i) {
                return new AccumulatorAverageBlockless();
            }
        });
        pSplitter.putDataInfo(meterP.getDataInfo());

        // collect pressure data before any insert/delete trials
        sim.integrator.addThermostatAction(new IAction() {
            public void actionPerformed() {
                pSplitter.putData(meterP.getData());
            }
        });

        final MCMoveIDBiasActionFasterer mcMoveBiasAction;
        final IntegratorListenerAction mcMoveBiasListener;

        if (params.doReweight) {
            mcMoveBiasAction = new MCMoveIDBiasActionFasterer(sim.integratorMC, sim.mcMoveID, numV,
                    mu, mcMoveOverlapMeter, numAtoms);
            mcMoveBiasAction.setNMaxReweight(numAtoms);
            mcMoveBiasAction.setDefaultDaDef(daDef);
            mcMoveBiasListener = new IntegratorListenerAction(mcMoveBiasAction, biasInterval);
            sim.integratorMC.getEventManager().addListener(mcMoveBiasListener);
        }
        else {
            mcMoveBiasAction = null;
            mcMoveBiasListener = null;
        }

        if (graphics) {
    	    
            MeterNMolecules meterN = new MeterNMolecules();
            meterN.setBox(sim.box);
            final AccumulatorHistogram nHistogram = new AccumulatorHistogram(new HistogramDiscrete(0.1));
            DataPumpListener nPump = new DataPumpListener(meterN, nHistogram, hybridInterval);
            sim.integrator.getEventManager().addListener(nPump);
            
            final String APP_NAME = "HSMD Vacancy";
        	final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);
            ColorScheme colorScheme = new ColorSchemeVacancy(sim.integrator, sim.neighborManager, sim.mcMoveID.getMaxDistance());

            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            
            DiameterHash dh = new DiameterHashVacancy(sim.integrator, sim.neighborManager, sim.mcMoveID.getMaxDistance());
            simGraphic.getDisplayBox(sim.box).setDiameterHash(dh);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
    
            final AccumulatorHistory nHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSourceCountTimeFasterer timeDataSource = new DataSourceCountTimeFasterer(sim.integrator);
            nHistory.setTimeDataSource(timeDataSource);
            nHistory.setPushInterval(1);
            DataFork nFork = new DataFork();
            nPump.setDataSink(nFork);
            nFork.addDataSink(nHistogram);
            nFork.addDataSink(nHistory);
            DisplayPlot nHistoryPlot = new DisplayPlot();
            nHistoryPlot.setLabel("N history");
            nHistoryPlot.setDoLegend(false);
            nHistory.setDataSink(nHistoryPlot.getDataSet().makeDataSink());
            simGraphic.add(nHistoryPlot);
            
            DisplayPlot nPlot = new DisplayPlot();
            nPlot.setLabel("N histogram");
            nPlot.getPlot().setYLog(true);
            nHistogram.setDataSink(nPlot.getDataSet().makeDataSink());
            nPlot.setLegend(new DataTag[]{nHistogram.getTag()}, "visited");
            simGraphic.add(nPlot);
            
            final DataSourceFEHistogramFasterer feHistogram = new DataSourceFEHistogramFasterer(mcMoveOverlapMeter, mu);
            DataPumpListener feHistPump = new DataPumpListener(feHistogram, nPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(feHistPump);
            nPlot.setLegend(new DataTag[]{feHistogram.getTag()}, "measured");
            
            final DataSourceImposedFEHistogramFasterer feHistogramImposed = new DataSourceImposedFEHistogramFasterer(mcMoveOverlapMeter, sim.mcMoveID, mu);
            DataPumpListener feHistImposedPump = new DataPumpListener(feHistogramImposed, nPlot.getDataSet().makeDataSink(), hybridInterval);
            nPlot.setLegend(new DataTag[]{feHistogramImposed.getTag()}, "imposed");
            sim.integrator.getEventManager().addListener(feHistImposedPump);

            if (!fluid) {
                final DataSourceImposedFEFitFasterer feHistogramFit = new DataSourceImposedFEFitFasterer(mcMoveOverlapMeter, sim.mcMoveID, mu);
                DataPumpListener feHistFitPump = new DataPumpListener(feHistogramFit, nPlot.getDataSet().makeDataSink(), hybridInterval);
                nPlot.setLegend(new DataTag[]{feHistogramFit.getTag()}, "fit");
                sim.integrator.getEventManager().addListener(feHistFitPump);
            }
            
            DisplayPlot fePlot = new DisplayPlot();
            final DataSourceFEFasterer dsfe = new DataSourceFEFasterer(mcMoveOverlapMeter);
            DataPumpListener fePump = new DataPumpListener(dsfe, fePlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(fePump);
            fePlot.setLegend(new DataTag[]{dsfe.getTag()}, "FE");

            if (!fluid) {
                final DataSourceFEFasterer dsfe2 = new DataSourceFEFasterer(mcMoveOverlapMeter);
                dsfe2.setSubtractComb(true);
                DataPumpListener fe2Pump = new DataPumpListener(dsfe2, fePlot.getDataSet().makeDataSink(), hybridInterval);
                sim.integrator.getEventManager().addListener(fe2Pump);
                fePlot.setLegend(new DataTag[]{dsfe2.getTag()}, "FEdef");
                
                final DataSourceFEFasterer dsfe3 = new DataSourceFEFasterer(mcMoveOverlapMeter);
                dsfe3.setSubtractComb(true);
                dsfe3.setAvgDef(true);
                DataPumpListener fe3Pump = new DataPumpListener(dsfe3, fePlot.getDataSet().makeDataSink(), hybridInterval);
                sim.integrator.getEventManager().addListener(fe3Pump);
                fePlot.setLegend(new DataTag[]{dsfe3.getTag()}, "FEdefAvg");
            }
            
            fePlot.setLabel("FE");
            simGraphic.add(fePlot);

            
            DisplayTextBox pDisplayBox = new DisplayTextBox();
            pDisplayBox.setLabel("Measured pressure");
            pDisplayBox.setPrecision(5);
            
            DisplayTextBox pDisplayBox2 = new DisplayTextBox();
            pDisplayBox2.setLabel("Pressure");
            pDisplayBox2.setPrecision(5);
            
            final DataSourceAvgPressureFasterer avgP = new DataSourceAvgPressureFasterer(pSplitter, mcMoveOverlapMeter, mu);
            DataPumpListener avgPPump = new DataPumpListener(avgP, pDisplayBox);
            sim.integrator.getEventManager().addListener(avgPPump);

            final DataSourceAvgPressure2Fasterer avgP2 = new DataSourceAvgPressure2Fasterer(pSplitter, mcMoveOverlapMeter, density, numAtoms/density);
            DataPumpListener avgPPump2 = new DataPumpListener(avgP2, pDisplayBox2);
            sim.integrator.getEventManager().addListener(avgPPump2);
            
            simGraphic.add(pDisplayBox);
            simGraphic.add(pDisplayBox2);

            DisplayPlot pmuPlot = new DisplayPlot();
            pmuPlot.setDoLegend(false);
            pmuPlot.setLabel("dP vs. mu");
            final DataSourcePMuFasterer pmu = new DataSourcePMuFasterer(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, false);
            DataPumpListener pmuPump = new DataPumpListener(pmu, pmuPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(pmuPump);
            
            simGraphic.add(pmuPlot);

            DisplayPlot muMuPlot = new DisplayPlot();
            muMuPlot.setDoLegend(false);
            muMuPlot.setLabel("mu vs. mu");
            final DataSourcePMuFasterer muMu = new DataSourcePMuFasterer(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, true);
            DataPumpListener muMuPump = new DataPumpListener(muMu, muMuPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(muMuPump);
            
            simGraphic.add(muMuPlot);

            DeviceBox muBox = new DeviceBox(sim.getController());
            muBox.setPrecision(8);
            muBox.setEditable(true);
            muBox.setLabel("mu");
            muBox.setModifier(new Modifier() {

                public double getValue() {
                    return sim.mcMoveID.getMu();
                }
                
                public void setValue(double newValue) {
                    sim.mcMoveID.setMu(newValue);
                    feHistogram.setMu(newValue);
                    feHistogramImposed.setMu(newValue);
                    avgP.setMu(newValue);
                    pmu.setMu(newValue);
                    if (mcMoveBiasAction != null) {
                        mcMoveBiasAction.setMu(newValue);
                    }
                }
                
                public String getLabel() {
                    return "your mom";
                }
                
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }
            });

            
            simGraphic.add(muBox);

            DataSourceMuRootFasterer dsmr = new DataSourceMuRootFasterer(mcMoveOverlapMeter, mu, pSplitter, mu, density, numAtoms/density);
            avgP2.setDataSourceMuRoot(dsmr);
            final AccumulatorHistory dsmrHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmrHistory.setTimeDataSource(timeDataSource);
            DataPumpListener dsmrPump = new DataPumpListener(dsmr, dsmrHistory, hybridInterval);
            sim.integrator.getEventManager().addListener(dsmrPump);
            DisplayPlot dsmrPlot = new DisplayPlot();
            dsmrPlot.setLabel("mu*");
            dsmrPlot.setLegend(new DataTag[]{dsmr.getTag()},"avg");
            dsmrHistory.setDataSink(dsmrPlot.getDataSet().makeDataSink());
            simGraphic.add(dsmrPlot);

            DataSourceMuRoot1Fasterer dsmr1 = null;
            final AccumulatorHistory dsmr1History;
            if (!fluid) {
                dsmr1 = new DataSourceMuRoot1Fasterer(mcMoveOverlapMeter, mu, pSplitter, Alat, density, numAtoms/density);
                dsmr1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
                dsmr1History.setTimeDataSource(timeDataSource);
                DataPumpListener dsmr1Pump = new DataPumpListener(dsmr1, dsmr1History, hybridInterval);
                sim.integrator.getEventManager().addListener(dsmr1Pump);
                dsmr1History.setDataSink(dsmrPlot.getDataSet().makeDataSink());
                dsmrPlot.setLegend(new DataTag[]{dsmr1.getTag()},"v1");
            }
            else {
                dsmr1History = null;
            }
            
            DataSourcePN dataSourcePN = new DataSourcePN(pSplitter, numAtoms);
            DisplayPlot plotPN = new DisplayPlot();
            DataPumpListener pumpPN = new DataPumpListener(dataSourcePN, plotPN.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(pumpPN);
            plotPN.setLabel("P vs. N");
            plotPN.setDoLegend(false);
            simGraphic.add(plotPN);
            
            DataSourceMuRootVacancyConcentration dsmrvc = dsmr.new DataSourceMuRootVacancyConcentration();
            final AccumulatorHistory dsmrvcHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmrvcHistory.setTimeDataSource(timeDataSource);
            DataPumpListener dsmrvcPump = new DataPumpListener(dsmrvc, dsmrvcHistory, hybridInterval);
            sim.integrator.getEventManager().addListener(dsmrvcPump);
            DisplayPlot vPlot = new DisplayPlot();
            vPlot.setLabel("vacancies");
            dsmrvcHistory.setDataSink(vPlot.getDataSet().makeDataSink());
            vPlot.setLegend(new DataTag[]{dsmrvc.getTag()}, "avg");
            simGraphic.add(vPlot);

            final AccumulatorHistory dsmrvc1History;
            if (!fluid) {
                DataSourceMuRoot1Fasterer.DataSourceMuRootVacancyConcentration dsmrvc1 = dsmr1.new DataSourceMuRootVacancyConcentration();
                dsmrvc1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
                dsmrvc1History.setTimeDataSource(timeDataSource);
                DataPumpListener dsmrvc1Pump = new DataPumpListener(dsmrvc1, dsmrvc1History, hybridInterval);
                sim.integrator.getEventManager().addListener(dsmrvc1Pump);
                dsmrvc1History.setDataSink(vPlot.getDataSet().makeDataSink());
                vPlot.setLegend(new DataTag[]{dsmrvc1.getTag()}, "v1");
            }
            else {
                dsmrvc1History = null;
            }

            final DeviceButton resetButton = new DeviceButton(sim.getController());
            IAction resetAction = new IAction() {
                boolean enabled = true;
                public void actionPerformed() {
                    if (enabled) {
                        mcMoveOverlapMeter.reset();
                        sim.integrator.resetStepCount();
                        meterP.getDataAsScalar();
                        meterP.reset();
                        nHistogram.reset();
                        nHistory.reset();
                        dsmrHistory.reset();
                        if (dsmr1History != null) {
                            dsmr1History.reset();
                            dsmrvc1History.reset();
                        }
                        dsmrvcHistory.reset();
                        for (int i=0; i<pSplitter.getNumDataSinks(); i++) {
                            AccumulatorAverageBlockless avg = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
                            if (avg != null) avg.reset();
                        }
                        sim.integratorMC.getEventManager().removeListener(mcMoveBiasListener);
                        enabled = false;
                        resetButton.setLabel("Reweight");
                    }
                    else {
                        sim.integratorMC.getEventManager().addListener(mcMoveBiasListener);
                        enabled = true;
                        resetButton.setLabel("Reset");
                    }
                }
            };
            resetButton.setAction(resetAction);
            resetButton.setLabel("Reset");
            simGraphic.getPanel().controlPanel.add(resetButton.graphic(),SimulationPanel.getVertGBC());

            return;
    	}

        final DataSourceFEFasterer dsfe3 = new DataSourceFEFasterer(mcMoveOverlapMeter);
        dsfe3.setSubtractComb(true);
        dsfe3.setAvgDef(true);


        // equilibrate off the lattice
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 40));

        IData dsfe3Data = dsfe3.getData();
        double daDefAvg = dsfe3Data.getValue(0);
        System.out.println("very initial daDef = "+daDefAvg);
        
        // throw away our data
        mcMoveOverlapMeter.reset();
        sim.integrator.resetStepCount();
        for (int i=0; i<pSplitter.getNumDataSinks(); i++) {
            AccumulatorAverageBlockless avg = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
            if (avg != null) avg.reset();
        }
        sim.integrator.resetStepCount();
        meterP.reset();

        if (params.doReweight) {
            // update weights
            mcMoveBiasAction.actionPerformed();
            // and stop adjusting them
            sim.integratorMC.getEventManager().removeListener(mcMoveBiasListener);
        }

        final long finalSteps = steps;
        sim.integrator.getEventManager().addListener(new IntegratorListener() {
            boolean reenabled = false;
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                if (!reenabled && sim.integrator.getStepCount() >= finalSteps/40) {
                    sim.integratorMC.getEventManager().addListener(mcMoveBiasListener);
                    reenabled = true;
                }
            }
            
            public void integratorInitialized(IntegratorEvent e) {}
        });

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));

        if (params.doReweight) {
            dsfe3Data = dsfe3.getData();
            daDefAvg = dsfe3Data.getValue(0);
            System.out.println("using daDef = "+daDefAvg);
            // update weights
            mcMoveBiasAction.actionPerformed();
            // and stop adjusting them (permanently this time)
            sim.integratorMC.getEventManager().removeListener(mcMoveBiasListener);
        }

        // and throw away data again
        mcMoveOverlapMeter.reset();
        for (int i=0; i<pSplitter.getNumDataSinks(); i++) {
            AccumulatorAverageBlockless avg = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
            if (avg != null) avg.reset();
        }
        sim.integrator.resetStepCount();
        meterP.reset();

        // take real data
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.currentTimeMillis();

        System.out.println();

        DataSourceMuRootFasterer dsmr = new DataSourceMuRootFasterer(mcMoveOverlapMeter, mu, pSplitter, mu, density, numAtoms/density);
        double muRoot = dsmr.getDataAsScalar();

        // histogram of # of atoms based on free energy differences
        final DataSourceFEHistogramFasterer feHistogram = new DataSourceFEHistogramFasterer(mcMoveOverlapMeter, muRoot);
        // actual free energy differences
        final DataSourceFEFasterer dsfe = new DataSourceFEFasterer(mcMoveOverlapMeter);
        final DataSourceFEFasterer dsfe2 = new DataSourceFEFasterer(mcMoveOverlapMeter);
        // subtract off combinatorial entropy
        dsfe2.setSubtractComb(true);
        IData nData = feHistogram.getIndependentData(0);
        IData fenData = feHistogram.getData();
        IData dsfeData = dsfe.getData();
        IData dsfe2Data = dsfe2.getData();
        dsfe3Data = dsfe3.getData();
        double[] nHistogram = mcMoveOverlapMeter.getHistogram();
        for (int i=0; i<nData.getLength(); i++) {
            int n = (int)Math.round(nData.getValue(i));
            double pAvg = Double.NaN;
            if (numAtoms-n < pSplitter.getNumDataSinks()) {
                AccumulatorAverageBlockless pAcc = (AccumulatorAverageBlockless)pSplitter.getDataSink(numAtoms-n);
                pAvg = pAcc == null ? Double.NaN : pAcc.getData().getValue(pAcc.AVERAGE.index);
            }
            if (Math.round(nData.getValue(i)) < numAtoms) {
                System.out.printf("%6d %20.15e %20.15e %20.15e %20.15e %20.15e\n", n, nHistogram[i], fenData.getValue(i), pAvg, dsfeData.getValue(i), dsfe2Data.getValue(i));
            }
            else {
                System.out.printf("%6d %20.15e %20.15e %20.15e\n", n, nHistogram[i], fenData.getValue(i), pAvg);
            }
        }
        System.out.println("\nfinal daDef: "+dsfe3Data.getValue(0));

        System.out.println("mu root: "+muRoot);
        double pRoot = dsmr.getLastPressure();
        double vRoot = dsmr.getLastVacancyConcentration();
        System.out.println("pressure root: "+pRoot);
        System.out.println("vacancy concentration root: "+vRoot);

        
        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }
    
    public static class HSMDVParams extends ParameterBase {
        public int numAtoms = 13500;
        public double density = 0.5;
        public double tStep = 4;
        public long steps = 40000;
        public int hybridInterval = -1;
        public boolean graphics = false;
        public double mu = Double.NaN;
        public double Alat = Double.NaN;
        public int numV = 5;
        public boolean doReweight = true;
        public double daDef = Double.NaN;
        public boolean fluid = false;
    }
}
