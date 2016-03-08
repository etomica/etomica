/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.awt.Color;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHash;
import etomica.box.Box;
import etomica.data.AccumulatorAverageBlockless;
import etomica.data.AccumulatorHistogram;
import etomica.data.AccumulatorHistory;
import etomica.data.DataDistributer;
import etomica.data.DataFork;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorScheme;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceButton;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveIDBiasAction;
import etomica.integrator.mcmove.MCMoveInsertDeleteLatticeVacancy;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.normalmode.DataSourceMuRoot.DataSourceMuRootVacancyConcentration;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.statmech.LennardJones;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.util.HistogramDiscrete;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.HistoryCollapsingDiscard;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimLJVacancy extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate ai;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public IBox box;
    public P2LennardJones p2LJ;
    public P2SoftSphericalTruncated potential;
    public MCMoveVolume mcMoveVolume;
    public MCMoveInsertDeleteLatticeVacancy mcMoveID;
    

    public SimLJVacancy(final int numAtoms, double temperature, double density, double rc, final int fixedN, final int maxDN, final double bmu) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        double L = Math.pow(4.0/density, 1.0/3.0);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        PrimitiveCubic primitive = new PrimitiveCubic(space, n*L);
        int[] nCells = new int[]{n,n,n};
        IBoundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        BasisBigCell basis = new BasisBigCell(space, basisFCC, nCells);
        
        box.setBoundary(boundary);

        
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setTemperature(temperature);
        MCMoveAtom move = new MCMoveAtom(random, potentialMaster, space);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(move);
        
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
        
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        if (rc > 0.49*box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("rc is too large");
        }
        double nbr1 = L/Math.sqrt(2);
        double y = 1.25*nbr1; //nbr1+(L-nbr1)*0.6+0.06;

        p2LJ = new P2LennardJones(space, 1, 1);
        potential = new P2SoftSphericalTruncated(space, p2LJ, rc);
        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType,leafType});

        integrator.setBox(box);
        
        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
        // we just needed this to initial coordinates, to keep track of lattice positions of atoms
        coordinateDefinition = null;
        
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        potentialMaster.getNbrCellManager(box).assignCellAll();
        
        mcMoveID = new MCMoveInsertDeleteLatticeVacancy(potentialMaster, random, space, integrator, y, fixedN, maxDN);

        double x = (nbr1-1)/4.0;
        mcMoveID.setMaxInsertDistance(x);
        mcMoveID.makeFccVectors(nbr1);
        mcMoveID.setMu(bmu);
        mcMoveID.setSpecies(species);
        integrator.getMoveManager().addMCMove(mcMoveID);
        integrator.getMoveManager().setFrequency(mcMoveID, 0.1);
        integrator.getMoveEventManager().addListener(mcMoveID);
    }

    public static void main(String[] args) {

        LJParams params = new LJParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.numAtoms = 500;
            params.steps = 1000000;
            params.density = 1.0;
            params.temperature = 0.7;
            params.numV = 5;

        }

        final int numAtoms = params.numAtoms;
        final double density = params.density;
        final double temperature = params.temperature;
        final double rc = params.rc;
        long steps = params.steps;
        boolean graphics = params.graphics;
        int maxDN = (params.numV+1)/2;
        int fixedN = params.numAtoms - maxDN;
        double mu = params.mu;
        double Alat = params.Alat;
        System.out.println("Running N="+params.numAtoms+" at rho="+density);
        if (Double.isNaN(mu)) {
            // this should get us pretty close
            Alat = Math.log(density)-1.0 + (LennardJones.aResidualFcc(temperature, density))/temperature;
            double zLat = LennardJones.ZFcc(temperature, density);
            System.out.println("lattice pressure "+zLat);
            mu = (Alat + zLat);
            System.out.println("beta mu "+mu);
        }
        else if (Double.isNaN(Alat)) {
            // this should get us pretty close
            Alat = Math.log(density)-1.0 + (LennardJones.aResidualFcc(temperature, density))/temperature;
        }
        double daDef = params.daDef;
        if (Double.isNaN(daDef)) {
            // use a correlation.  this should be with ~1
            daDef = 28 - 20*temperature;
            System.out.println("estimated defect free energy: "+daDef);
        }
        else {
            System.out.println("using defect free energy: "+daDef);
        }
        System.out.println("Running LJ MC");
        System.out.println("N: "+numAtoms);
        System.out.println("T: "+temperature);
        System.out.println("density: "+density);
        System.out.println("steps: "+steps);
   

        final SimLJVacancy sim = new SimLJVacancy(numAtoms, temperature, density, rc, fixedN, maxDN, mu);


        final int biasInterval = 10*numAtoms;

        final MCMoveOverlapListener mcMoveOverlapMeter = new MCMoveOverlapListener(sim.mcMoveID, 11, daDef, numAtoms, 6);
        mcMoveOverlapMeter.setTemperature(temperature);
        sim.integrator.getMoveEventManager().addListener(mcMoveOverlapMeter);

        final MeterPressure meterP = new MeterPressure(sim.getSpace());
        meterP.setIntegrator(sim.integrator);
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
        sim.integrator.getEventManager().addListener(new IIntegratorListener() {
            long countDown = numAtoms, interval = numAtoms;
            public void integratorInitialized(IIntegratorEvent e) {}

            public void integratorStepStarted(IIntegratorEvent e) {}

            public void integratorStepFinished(IIntegratorEvent e) {
                countDown--;
                if (countDown==0) {
                    pSplitter.putData(meterP.getData());
                    countDown=interval;
                }
            }
        });

        final MCMoveIDBiasAction mcMoveBiasAction;
        final IntegratorListenerAction mcMoveBiasListener;

        if (params.doReweight) {
            mcMoveBiasAction = new MCMoveIDBiasAction(sim.integrator, sim.mcMoveID, maxDN,
                    fixedN, mu, mcMoveOverlapMeter, numAtoms, 1);
            mcMoveBiasAction.setNMaxReweight(numAtoms);
            mcMoveBiasAction.setDefaultDaDef(daDef);
            mcMoveBiasListener = new IntegratorListenerAction(mcMoveBiasAction, biasInterval);
            sim.integrator.getEventManager().addListener(mcMoveBiasListener);
        }
        else {
            mcMoveBiasAction = null;
            mcMoveBiasListener = null;
        }

        if (graphics) {
    	    
            MeterNMolecules meterN = new MeterNMolecules();
            meterN.setBox(sim.box);
            final AccumulatorHistogram nHistogram = new AccumulatorHistogram(new HistogramDiscrete(0.1));
            DataPumpListener nPump = new DataPumpListener(meterN, nHistogram, numAtoms);
            sim.integrator.getEventManager().addListener(nPump);
            
            final String APP_NAME = "LJ Vacancy";
        	final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, numAtoms, sim.space, sim.getController());
        	final int nc = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            ColorScheme colorScheme = new ColorScheme() {
                public Color getAtomColor(IAtom a) {
                    float x = (float)a.getPosition().getX(0);
                    float L = (float)sim.box.getBoundary().getBoxSize().getX(0);
                    x = nc*2*(x/L+0.5f);
                    if (nc%2==0) x+=0.5f;
                    while (x>3) x-=3;
                    float r = 0f, g = 0f, b = 0f;
                    if (x < 1) {
                        r = 1-x;
                        if (r>1) r=1;
                        g = x;
                        if (g<0) g=0;
                    }
                    else if (x < 2) {
                        g = 2-x;
                        b = x-1;
                    }
                    else {
                        b = 3-x;
                        if (b<0) b=0;
                        r = x-2;
                        if (r>1) r=1;
                    }
                    return new Color(r, g, b);
                }
            };
            colorScheme = new ColorScheme() {
                IVectorMutable dr = sim.space.makeVector();
                double rMax = sim.mcMoveID.getMaxDistance();
                double rc2 = rMax*rMax;
                int nmax = 12;
                Api1ACell iter = new Api1ACell(3, rMax, sim.potentialMaster.getCellAgentManager());
                public Color getAtomColor(IAtom a) {
                    if (!sim.integrator.getEventManager().firingEvent() && !sim.ai.isPaused()) return new Color(1.0f, 1.0f, 1.0f);

                    IVector pi = a.getPosition();
                    iter.setBox(sim.box);
                    iter.setDirection(null);
                    iter.setTarget(a);
                    iter.reset();
                    IBoundary boundary = sim.box.getBoundary();
                    int n = 0;
                    for (IAtomList pair = iter.next(); pair!=null; pair=iter.next()) {
                        IAtom nbrj = pair.getAtom(0);
                        if (nbrj==a) nbrj=pair.getAtom(1);
                        dr.Ev1Mv2(pi, nbrj.getPosition());
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();
                        if (r2 < rc2) {
                            n++;
                        }
                    }
                    if (n < nmax) return Color.RED;
                    if (n == nmax) return Color.BLUE;
                    return Color.BLUE;
                }
            };

            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            
            DiameterHash dh = new DiameterHash() {
                IVectorMutable dr = sim.space.makeVector();
                double rMax = sim.mcMoveID.getMaxDistance();
                double rc2 = rc*rc;
                int nmax = 12;
                Api1ACell iter = new Api1ACell(3, rMax, sim.potentialMaster.getCellAgentManager());
                public double getDiameter(IAtom a) {
                    IVector pi = a.getPosition();
                    IBoundary boundary = sim.box.getBoundary();
                    iter.setBox(sim.box);
                    iter.setDirection(null);
                    iter.setTarget(a);
                    try {
                        iter.reset();
                    }
                    catch (NullPointerException ex) {
                        // during addition
                        return 0;
                    }
                    int n = 0;
                    for (IAtomList pair = iter.next(); pair!=null; pair=iter.next()) {
                        IAtom nbrj = pair.getAtom(0);
                        if (nbrj==a) nbrj=pair.getAtom(1);
                        dr.Ev1Mv2(pi, nbrj.getPosition());
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();
                        if (r2 < rc2) {
                            n++;
                        }
                    }
                    if (n < nmax) return 1.0;
                    if (n == nmax) return 0.1;
                    return 0.1;
                }
            };
//            simGraphic.getDisplayBox(sim.box).setDiameterHash(dh);

            
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
    
            final AccumulatorHistory nHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSourceCountSteps timeDataSource = new DataSourceCountSteps(sim.integrator);
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
            
            final DataSourceFEHistogram feHistogram = new DataSourceFEHistogram(mcMoveOverlapMeter, mu);
            DataPumpListener feHistPump = new DataPumpListener(feHistogram, nPlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(feHistPump);
            nPlot.setLegend(new DataTag[]{feHistogram.getTag()}, "measured");
            
            final DataSourceImposedFEHistogram feHistogramImposed = new DataSourceImposedFEHistogram(mcMoveOverlapMeter, sim.mcMoveID, mu);
            DataPumpListener feHistImposedPump = new DataPumpListener(feHistogramImposed, nPlot.getDataSet().makeDataSink(), numAtoms);
            nPlot.setLegend(new DataTag[]{feHistogramImposed.getTag()}, "imposed");
            sim.integrator.getEventManager().addListener(feHistImposedPump);

            final DataSourceImposedFEFit feHistogramFit = new DataSourceImposedFEFit(mcMoveOverlapMeter, sim.mcMoveID, mu);
            DataPumpListener feHistFitPump = new DataPumpListener(feHistogramFit, nPlot.getDataSet().makeDataSink(), numAtoms);
            nPlot.setLegend(new DataTag[]{feHistogramFit.getTag()}, "fit");
            sim.integrator.getEventManager().addListener(feHistFitPump);
            
            DisplayPlot fePlot = new DisplayPlot();
            final DataSourceFE dsfe = new DataSourceFE(mcMoveOverlapMeter);
            DataPumpListener fePump = new DataPumpListener(dsfe, fePlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(fePump);
            fePlot.setLegend(new DataTag[]{dsfe.getTag()}, "FE");

            final DataSourceFE dsfe2 = new DataSourceFE(mcMoveOverlapMeter);
            dsfe2.setSubtractComb(true);
            DataPumpListener fe2Pump = new DataPumpListener(dsfe2, fePlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(fe2Pump);
            fePlot.setLegend(new DataTag[]{dsfe2.getTag()}, "FEdef");
            
            final DataSourceFE dsfe3 = new DataSourceFE(mcMoveOverlapMeter);
            dsfe3.setSubtractComb(true);
            dsfe3.setAvgDef(true);
            DataPumpListener fe3Pump = new DataPumpListener(dsfe3, fePlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(fe3Pump);
            fePlot.setLegend(new DataTag[]{dsfe3.getTag()}, "FEdefAvg");
            
            fePlot.setLabel("FE");
            simGraphic.add(fePlot);

            
            DisplayTextBox pDisplayBox = new DisplayTextBox();
            pDisplayBox.setLabel("Measured pressure");
            pDisplayBox.setPrecision(5);
            
            DisplayTextBox pDisplayBox2 = new DisplayTextBox();
            pDisplayBox2.setLabel("Pressure");
            pDisplayBox2.setPrecision(5);
            
            final DataSourceAvgPressure avgP = new DataSourceAvgPressure(pSplitter, mcMoveOverlapMeter, mu);
            DataPumpListener avgPPump = new DataPumpListener(avgP, pDisplayBox, numAtoms);
            sim.integrator.getEventManager().addListener(avgPPump);

            final DataSourceAvgPressure2 avgP2 = new DataSourceAvgPressure2(mcMoveOverlapMeter, mu, Alat, density, numAtoms/density);
            DataPumpListener avgPPump2 = new DataPumpListener(avgP2, pDisplayBox2, numAtoms);
            sim.integrator.getEventManager().addListener(avgPPump2);
            
            simGraphic.add(pDisplayBox);
            simGraphic.add(pDisplayBox2);

            DisplayPlot pmuPlot = new DisplayPlot();
            pmuPlot.setDoLegend(false);
            pmuPlot.setLabel("dP vs. mu");
            final DataSourcePMu pmu = new DataSourcePMu(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, false);
            DataPumpListener pmuPump = new DataPumpListener(pmu, pmuPlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(pmuPump);
            
            simGraphic.add(pmuPlot);

            DisplayPlot muMuPlot = new DisplayPlot();
            muMuPlot.setDoLegend(false);
            muMuPlot.setLabel("mu vs. mu");
            final DataSourcePMu muMu = new DataSourcePMu(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, true);
            DataPumpListener muMuPump = new DataPumpListener(muMu, muMuPlot.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(muMuPump);
            
            simGraphic.add(muMuPlot);

            DeviceBox muBox = new DeviceBox();
            muBox.setPrecision(8);
            muBox.setEditable(true);
            muBox.setLabel("mu");
            muBox.setModifier(new Modifier() {
                
                public void setValue(double newValue) {
                    sim.mcMoveID.setMu(newValue);
                    feHistogram.setMu(newValue);
                    feHistogramImposed.setMu(newValue);
                    avgP.setMu(newValue);
                    avgP2.setBetaMu(newValue);
                    pmu.setMu(newValue);
                    mcMoveBiasAction.setMu(newValue);
                }
                
                public double getValue() {
                    return sim.mcMoveID.getMu();
                }
                
                public String getLabel() {
                    return "your mom";
                }
                
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }
            });

            
            simGraphic.add(muBox);

            DataSourceMuRoot dsmr = new DataSourceMuRoot(mcMoveOverlapMeter, mu, pSplitter, Alat, density, numAtoms/density);
            final AccumulatorHistory dsmrHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmrHistory.setTimeDataSource(timeDataSource);
            DataPumpListener dsmrPump = new DataPumpListener(dsmr, dsmrHistory, numAtoms);
            sim.integrator.getEventManager().addListener(dsmrPump);
            DisplayPlot dsmrPlot = new DisplayPlot();
            dsmrPlot.setLabel("mu*");
            dsmrPlot.setLegend(new DataTag[]{dsmr.getTag()},"avg");
            dsmrHistory.setDataSink(dsmrPlot.getDataSet().makeDataSink());
            simGraphic.add(dsmrPlot);

            DataSourceMuRoot1 dsmr1 = null;
            final AccumulatorHistory dsmr1History;
            dsmr1 = new DataSourceMuRoot1(mcMoveOverlapMeter, mu, pSplitter, Alat, density, numAtoms/density);
            dsmr1.setMinMu(10);
            dsmr1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmr1History.setTimeDataSource(timeDataSource);
            DataPumpListener dsmr1Pump = new DataPumpListener(dsmr1, dsmr1History, numAtoms);
            sim.integrator.getEventManager().addListener(dsmr1Pump);
            dsmr1History.setDataSink(dsmrPlot.getDataSet().makeDataSink());
            dsmrPlot.setLegend(new DataTag[]{dsmr1.getTag()},"v1");
            
            DataSourcePN dataSourcePN = new DataSourcePN(pSplitter, numAtoms);
            DisplayPlot plotPN = new DisplayPlot();
            DataPumpListener pumpPN = new DataPumpListener(dataSourcePN, plotPN.getDataSet().makeDataSink(), numAtoms);
            sim.integrator.getEventManager().addListener(pumpPN);
            plotPN.setLabel("P vs. N");
            plotPN.setDoLegend(false);
            simGraphic.add(plotPN);
            
            DataSourceMuRootVacancyConcentration dsmrvc = dsmr.new DataSourceMuRootVacancyConcentration();
            final AccumulatorHistory dsmrvcHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmrvcHistory.setTimeDataSource(timeDataSource);
            DataPumpListener dsmrvcPump = new DataPumpListener(dsmrvc, dsmrvcHistory, numAtoms);
            sim.integrator.getEventManager().addListener(dsmrvcPump);
            DisplayPlot vPlot = new DisplayPlot();
            vPlot.setLabel("vacancies");
            dsmrvcHistory.setDataSink(vPlot.getDataSet().makeDataSink());
            vPlot.setLegend(new DataTag[]{dsmrvc.getTag()}, "avg");
            simGraphic.add(vPlot);

            final AccumulatorHistory dsmrvc1History;
            DataSourceMuRoot1.DataSourceMuRootVacancyConcentration dsmrvc1 = dsmr1.new DataSourceMuRootVacancyConcentration();
            dsmrvc1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmrvc1History.setTimeDataSource(timeDataSource);
            DataPumpListener dsmrvc1Pump = new DataPumpListener(dsmrvc1, dsmrvc1History, numAtoms);
            sim.integrator.getEventManager().addListener(dsmrvc1Pump);
            dsmrvc1History.setDataSink(vPlot.getDataSet().makeDataSink());
            vPlot.setLegend(new DataTag[]{dsmrvc1.getTag()}, "v1");

            final DeviceButton resetButton = new DeviceButton(sim.getController());
            IAction resetAction = new IAction() {
                boolean enabled = true;
                public void actionPerformed() {
                    if (enabled) {
                        mcMoveOverlapMeter.reset();
                        sim.integrator.resetStepCount();
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
                        sim.integrator.getEventManager().removeListener(mcMoveBiasListener);
                        enabled = false;
                        resetButton.setLabel("Reweight");
                    }
                    else {
                        sim.integrator.getEventManager().addListener(mcMoveBiasListener);
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

        final DataSourceFE dsfe3 = new DataSourceFE(mcMoveOverlapMeter);
        dsfe3.setSubtractComb(true);
        dsfe3.setAvgDef(true);


        // equilibrate off the lattice
        sim.ai.setMaxSteps(steps/40);
        sim.ai.actionPerformed();

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

        if (params.doReweight) {
            // update weights
            mcMoveBiasAction.actionPerformed();
            // and stop adjusting them
            sim.integrator.getEventManager().removeListener(mcMoveBiasListener);
        }

        final long finalSteps = steps;
        sim.integrator.getEventManager().addListener(new IIntegratorListener() {
            boolean reenabled = false;
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
                if (!reenabled && sim.integrator.getStepCount() >= finalSteps/40) {
                    sim.integrator.getEventManager().addListener(mcMoveBiasListener);
                    reenabled = true;
                }
            }
            
            public void integratorInitialized(IIntegratorEvent e) {}
        });

        sim.ai.setMaxSteps(steps/10);
        sim.ai.actionPerformed();

        if (params.doReweight) {
            dsfe3Data = dsfe3.getData();
            daDefAvg = dsfe3Data.getValue(0);
            System.out.println("using daDef = "+daDefAvg);
            // update weights
            mcMoveBiasAction.actionPerformed();
            // and stop adjusting them (permanently this time)
            sim.integrator.getEventManager().removeListener(mcMoveBiasListener);
        }

        // and throw away data again
        mcMoveOverlapMeter.reset();
        for (int i=0; i<pSplitter.getNumDataSinks(); i++) {
            AccumulatorAverageBlockless avg = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
            if (avg != null) avg.reset();
        }
        sim.integrator.resetStepCount();

        
        // take real data
        long t1 = System.currentTimeMillis();
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        System.out.println();

        DataSourceMuRoot dsmr = new DataSourceMuRoot(mcMoveOverlapMeter, mu, pSplitter, Alat, density, numAtoms/density);
        double muRoot = dsmr.getDataAsScalar();

        // histogram of # of atoms based on free energy differences
        final DataSourceFEHistogram feHistogram = new DataSourceFEHistogram(mcMoveOverlapMeter, muRoot);
        // actual free energy differences
        final DataSourceFE dsfe = new DataSourceFE(mcMoveOverlapMeter);
        final DataSourceFE dsfe2 = new DataSourceFE(mcMoveOverlapMeter);
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
                System.out.println(String.format("%6d %20.15e %20.15e %20.15e %20.15e %20.15e", n, nHistogram[i], fenData.getValue(i), pAvg, dsfeData.getValue(i), dsfe2Data.getValue(i)));
            }
            else {
                System.out.println(String.format("%6d %20.15e %20.15e %20.15e", n, nHistogram[i], fenData.getValue(i), pAvg));
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
    
    public static class LJParams extends ParameterBase {
        public int numAtoms = 500;
        public double density = 1.0;
        public double temperature = 0.8;
        public long steps = 1000000;
        public boolean graphics = false;
        public double mu = Double.NaN;
        public double Alat = Double.NaN;
        public int numV = 5;
        public boolean doReweight = true;
        public double daDef = Double.NaN;
        public double rc = 3;
    }
}
