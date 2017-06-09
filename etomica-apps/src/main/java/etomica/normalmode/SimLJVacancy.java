/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.space.Boundary;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.normalmode.DataSourceMuRoot.DataSourceMuRootVacancyConcentration;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.statmech.LennardJones;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimLJVacancy extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate ai;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public Potential2SoftSpherical p2LJ;
    public P2SoftSphericalTruncated potential;
    public MCMoveVolume mcMoveVolume;
    public MCMoveInsertDeleteLatticeVacancy mcMoveID;
    public double uShift;

    public SimLJVacancy(final int numAtoms, double temperature, double density, double rc, final int numV, final double bmu, final boolean ss, boolean shifted) {
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
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        BasisBigCell basis = new BasisBigCell(space, basisFCC, nCells);
        
        box.setBoundary(boundary);

        
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setTemperature(temperature);
        MCMoveAtom move = new MCMoveAtom(random, potentialMaster, space);
        move.setStepSize(0.2);
//        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(move);
        
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
        
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        integrator.setBox(box);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
        // we just needed this to initial coordinates, to keep track of lattice positions of atoms
        coordinateDefinition = null;

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        if (rc > 0.49*box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("rc is too large");
        }
        p2LJ = ss ? new P2SoftSphere(space, 1, 4, 12) : new P2LennardJones(space, 1, 1);
        potential = new P2SoftSphericalTruncated(space, p2LJ, rc);
        potential.setMakeLrc(false);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        potentialMaster.getNbrCellManager(box).assignCellAll();
        if (shifted) {
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
            meterPE.setBox(box);
            meterPE.setIncludeLrc(false);
            double uUnshifted = meterPE.getDataAsScalar();
            
            potentialMaster.removePotential(potential);
            potential = new P2SoftSphericalTruncatedShifted(space, p2LJ, rc);
            potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});
            double uShifted = meterPE.getDataAsScalar();
            uShift = (uUnshifted - uShifted)/numAtoms;
            Potential0Lrc pShift = new Potential0Lrc(space, new AtomType[]{leafType, leafType}, potential) {
                public double virial(IAtomList atoms) {
                    return 0;
                }
                
                public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
                    return null;
                }
                
                public Vector[] gradient(IAtomList atoms) {
                    return null;
                }
                
                public double energy(IAtomList atoms) {
                    if (divisor==0) return 0;
                    if (divisor==1) {
                        // the whole system
                        return numAtoms*uShift;
                    }
                    // single atom
                    return 2*uShift;
                }
            };
            potentialMaster.lrcMaster().addPotential(pShift);
        }

        double nbr1 = L/Math.sqrt(2);
        double y = 1.25*nbr1;

        mcMoveID = new MCMoveInsertDeleteLatticeVacancy(potentialMaster, random, space, integrator, y, numAtoms, numV);

        mcMoveID.setMaxInsertDistance(nbr1/40);
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
            params.density = 1.5;
            params.numV = 1;
            params.temperature = 1;
        }

        final int numAtoms = params.numAtoms;
        final double density = params.density;
        final double temperature = params.temperature;
        final double rc = params.rc*Math.pow(density, -1.0/3.0);
        long steps = params.steps;
        boolean graphics = params.graphics;
        int numV = params.numV;
        double mu = params.mu;
        boolean fixedDaDef = params.fixedDaDef;
        boolean ss = params.ss;
        boolean shifted = params.shifted;
        System.out.println("Running N="+params.numAtoms+" at rho="+density);
        if (Double.isNaN(mu)) {
            // this should get us pretty close
            double Alat = Math.log(density)-1.0 + (LennardJones.aResidualFcc(temperature, density))/temperature;
            double zLat = LennardJones.ZFcc(temperature, density);
            mu = (Alat + zLat);
        }
        System.out.println("beta mu "+mu);
        double daDef = params.daDef;
        if (Double.isNaN(daDef)) {
            // use a correlation.  this should be with ~1
            double Y = temperature/Math.pow(density, 4)/4;
            double v2 = 1/(density*density);
            if (ss) v2 = 0;
            daDef = (-3.742    - 9.925*Y     - 2.081*Y*Y
                     +8.209*v2 + 1.929*Y*v2
                     -0.3479*v2*v2)/Y;
            System.out.println("estimated defect free energy: "+daDef);
        }
        else {
            System.out.println("using defect free energy: "+daDef);
        }
        System.out.println("Running "+(ss?"SS":"LJ")+" MC");
        System.out.println("N: "+numAtoms);
        System.out.println("T: "+temperature);
        System.out.println("density: "+density);
        System.out.println("steps: "+steps);
   

        final SimLJVacancy sim = new SimLJVacancy(numAtoms, temperature, density, rc, numV, mu, ss, shifted);
        if (shifted) {
            System.out.println("Shifting the potential (and then correcting by "+sim.uShift+")");
        }

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
        sim.integrator.getEventManager().addListener(new IntegratorListener() {
            long countDown = numAtoms, interval = numAtoms;
            public void integratorInitialized(IntegratorEvent e) {}

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                countDown--;
                if (countDown==0) {
                    // everything really wants beta P
                    IData pData = meterP.getData();
                    pData.TE(1.0/temperature);
                    pSplitter.putData(pData);
                    countDown=interval;
                }
            }
        });

        final MCMoveIDBiasAction mcMoveBiasAction;
        final IntegratorListenerAction mcMoveBiasListener;

        if (params.doReweight) {
            mcMoveBiasAction = new MCMoveIDBiasAction(sim.integrator, sim.mcMoveID, numV,
                    mu, mcMoveOverlapMeter, numAtoms);
            mcMoveBiasAction.setNMaxReweight(numAtoms);
            mcMoveBiasAction.setDefaultDaDef(daDef);
            if (fixedDaDef) {
                mcMoveBiasAction.setFixedDaDef(daDef);
            }
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
                Vector dr = sim.space.makeVector();
                double rMax = sim.mcMoveID.getMaxDistance();
                double rc2 = rMax*rMax;
                int nmax = 12;
                Api1ACell iter = new Api1ACell(3, rMax, sim.potentialMaster.getCellAgentManager());
                public Color getAtomColor(IAtom a) {
                    if (!sim.integrator.getEventManager().firingEvent() && !sim.ai.isPaused()) return new Color(1.0f, 1.0f, 1.0f);

                    Vector pi = a.getPosition();
                    iter.setBox(sim.box);
                    iter.setDirection(null);
                    iter.setTarget(a);
                    iter.reset();
                    Boundary boundary = sim.box.getBoundary();
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

            final DataSourceAvgPressure2 avgP2 = new DataSourceAvgPressure2(pSplitter, mcMoveOverlapMeter, density, numAtoms/density);
            DataPumpListener avgPPump2 = new DataPumpListener(avgP2, pDisplayBox2, numAtoms);
            sim.integrator.getEventManager().addListener(avgPPump2);
            
            simGraphic.add(pDisplayBox);
            simGraphic.add(pDisplayBox2);

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

            DataSourceMuRoot dsmr = new DataSourceMuRoot(mcMoveOverlapMeter, mu, pSplitter, mu, density, numAtoms/density);
            avgP2.setDataSourceMuRoot(dsmr);
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
            dsmr1 = new DataSourceMuRoot1(mcMoveOverlapMeter, mu, pSplitter, Double.NaN, density, numAtoms/density);
            dsmr1.setMinMu(10);
            dsmr1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmr1History.setTimeDataSource(timeDataSource);
            DataPumpListener dsmr1Pump = new DataPumpListener(dsmr1, dsmr1History, numAtoms);
//            sim.integrator.getEventManager().addListener(dsmr1Pump);
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
//            sim.integrator.getEventManager().addListener(dsmrvc1Pump);
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


        double daDefAvg = 0;
        IData dsfe3Data = null;
        if (!fixedDaDef) {
            // equilibrate off the lattice
            sim.ai.setMaxSteps(steps/40);
            sim.ai.actionPerformed();
    
            dsfe3Data = dsfe3.getData();
            if (dsfe3Data.getLength() == 0) {
                throw new RuntimeException("no data after initial equilibration");
            }
            daDefAvg = dsfe3Data.getValue(0);
            System.out.println("very initial daDef = "+daDefAvg);
            
            // throw away our data
            mcMoveOverlapMeter.reset();
            sim.integrator.resetStepCount();
            for (int i=0; i<pSplitter.getNumDataSinks(); i++) {
                AccumulatorAverageBlockless avg = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
                if (avg != null) avg.reset();
            }
            sim.integrator.resetStepCount();
        }

        if (fixedDaDef || params.doReweight) {
            // update weights
            mcMoveBiasAction.actionPerformed();
            // and stop adjusting them
            sim.integrator.getEventManager().removeListener(mcMoveBiasListener);
        }

        final long finalSteps = steps;
        if (!fixedDaDef) {
            // wait until 1/4 of the way through 2nd stage initialization, then start readjusting weights again
            sim.integrator.getEventManager().addListener(new IntegratorListener() {
                boolean reenabled = false;
                public void integratorStepStarted(IntegratorEvent e) {}

                public void integratorStepFinished(IntegratorEvent e) {
                    if (!reenabled && sim.integrator.getStepCount() >= finalSteps/40) {
                        sim.integrator.getEventManager().addListener(mcMoveBiasListener);
                        reenabled = true;
                    }
                }

                public void integratorInitialized(IntegratorEvent e) {}
            });
        }

        sim.ai.setMaxSteps(steps/10);
        sim.ai.actionPerformed();
        System.out.println("equilibration finished");

        if (!fixedDaDef && params.doReweight) {
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
        sim.integrator.getMoveManager().setEquilibrating(false);

        
        // take real data
        long t1 = System.currentTimeMillis();
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        System.out.println();

        DataSourceMuRoot dsmr = new DataSourceMuRoot(mcMoveOverlapMeter, mu, pSplitter, mu, density, numAtoms/density);
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
                pAvg = pAcc == null ? Double.NaN : pAcc.getData().getValue(AccumulatorAverageBlockless.AVERAGE.index);
            }
            if (dsfeData.getLength() > i) {
                System.out.println(String.format("%6d %20.15e %20.15e %20.15e %20.15e %20.15e", n, nHistogram[i], fenData.getValue(i), pAvg, dsfeData.getValue(i), dsfe2Data.getValue(i)));
            }
            else {
                System.out.println(String.format("%6d %20.15e %20.15e %20.15e", n, nHistogram[i], fenData.getValue(i), pAvg));
            }
        }

        
        double pAccept = sim.mcMoveID.getTracker().acceptanceProbability();
        System.out.println("\nInsert/delete acceptance "+pAccept);

        System.out.println("final daDef: "+dsfe3Data.getValue(0));

        System.out.println("bmu: "+muRoot);
        double pRoot = dsmr.getLastPressure();
        double vRoot = dsmr.getLastVacancyConcentration();
        System.out.println("pressure: "+pRoot);
        double deltaP = dsmr.getLastDPressure();
        System.out.println("delta P: "+deltaP);
        double deltaMu = dsmr.getLastDMu();
        System.out.println("delta bmu: "+deltaMu);
        System.out.println("vacancy fraction: "+vRoot);
        double lnTot = dsmr.getLastLnTot();
        double deltaA = -(1-vRoot)*lnTot/numAtoms - vRoot/density*deltaP/temperature;
        // free energy change per lattice site
        System.out.println("delta betaA: "+deltaA);
        
        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }
    
    public static class Applet extends javax.swing.JApplet {
    	public void init() {
            LJParams params = new LJParams();
            params.graphics = true;
            params.numAtoms = 864;
            params.steps = 1000000;
            params.density = 1.;
            params.numV = 4;
            params.temperature = 1.08;
            params.mu = 1.741248580731481;
            params.daDef = 6.39953903744669;
            params.rc = 4;

            final int numAtoms = params.numAtoms;
            final double density = params.density;
            final double temperature = params.temperature;
            final double rc = params.rc*Math.pow(density, -1.0/3.0);
            long steps = params.steps;
            int numV = params.numV;
            double mu = params.mu;
            boolean fixedDaDef = params.fixedDaDef;
            boolean ss = params.ss;
            boolean shifted = params.shifted;
            System.out.println("Running N="+params.numAtoms+" at rho="+density);
            if (Double.isNaN(mu)) {
                // this should get us pretty close
                double Alat = Math.log(density)-1.0 + (LennardJones.aResidualFcc(temperature, density))/temperature;
                double zLat = LennardJones.ZFcc(temperature, density);
                mu = (Alat + zLat);
            }
            System.out.println("beta mu "+mu);
            double daDef = params.daDef;
            if (Double.isNaN(daDef)) {
                // use a correlation.  this should be with ~1
                double Y = temperature/Math.pow(density, 4)/4;
                double v2 = 1/(density*density);
                if (ss) v2 = 0;
                daDef = (-3.742    - 9.925*Y     - 2.081*Y*Y
                         +8.209*v2 + 1.929*Y*v2
                         -0.3479*v2*v2)/Y;
                System.out.println("estimated defect free energy: "+daDef);
            }
            else {
                System.out.println("using defect free energy: "+daDef);
            }
            System.out.println("Running "+(ss?"SS":"LJ")+" MC");
            System.out.println("N: "+numAtoms);
            System.out.println("T: "+temperature);
            System.out.println("density: "+density);
            System.out.println("steps: "+steps);
       

            final SimLJVacancy sim = new SimLJVacancy(numAtoms, temperature, density, rc, numV, mu, ss, shifted);
            if (shifted) {
                System.out.println("Shifting the potential (and then correcting by "+sim.uShift+")");
            }

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
            sim.integrator.getEventManager().addListener(new IntegratorListener() {
                long countDown = numAtoms, interval = numAtoms;
                public void integratorInitialized(IntegratorEvent e) {}

                public void integratorStepStarted(IntegratorEvent e) {}

                public void integratorStepFinished(IntegratorEvent e) {
                    countDown--;
                    if (countDown==0) {
                        // everything really wants beta P
                        IData pData = meterP.getData();
                        pData.TE(1.0/temperature);
                        pSplitter.putData(pData);
                        countDown=interval;
                    }
                }
            });

            final MCMoveIDBiasAction mcMoveBiasAction;
            final IntegratorListenerAction mcMoveBiasListener;

            if (params.doReweight) {
                mcMoveBiasAction = new MCMoveIDBiasAction(sim.integrator, sim.mcMoveID, numV,
                        mu, mcMoveOverlapMeter, numAtoms);
                mcMoveBiasAction.setNMaxReweight(numAtoms);
                mcMoveBiasAction.setDefaultDaDef(daDef);
                if (fixedDaDef) {
                    mcMoveBiasAction.setFixedDaDef(daDef);
                }
                mcMoveBiasListener = new IntegratorListenerAction(mcMoveBiasAction, biasInterval);
                sim.integrator.getEventManager().addListener(mcMoveBiasListener);
            }
            else {
                mcMoveBiasAction = null;
                mcMoveBiasListener = null;
            }

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
                Vector dr = sim.space.makeVector();
                double rMax = sim.mcMoveID.getMaxDistance();
                double rc2 = rMax*rMax;
                int nmax = 12;
                Api1ACell iter = new Api1ACell(3, rMax, sim.potentialMaster.getCellAgentManager());
                public Color getAtomColor(IAtom a) {
                    if (!sim.integrator.getEventManager().firingEvent() && !sim.ai.isPaused()) return new Color(1.0f, 1.0f, 1.0f);

                    Vector pi = a.getPosition();
                    iter.setBox(sim.box);
                    iter.setDirection(null);
                    iter.setTarget(a);
                    iter.reset();
                    Boundary boundary = sim.box.getBoundary();
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

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            
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

            final DataSourceAvgPressure2 avgP2 = new DataSourceAvgPressure2(pSplitter, mcMoveOverlapMeter, density, numAtoms/density);
            DataPumpListener avgPPump2 = new DataPumpListener(avgP2, pDisplayBox2, numAtoms);
            sim.integrator.getEventManager().addListener(avgPPump2);
            
            simGraphic.add(pDisplayBox);
            simGraphic.add(pDisplayBox2);

            DeviceBox muBox = new DeviceBox();
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
                    mcMoveBiasAction.setMu(newValue);
                }
                
                public String getLabel() {
                    return "your mom";
                }
                
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }
            });

            
            simGraphic.add(muBox);

            DataSourceMuRoot dsmr = new DataSourceMuRoot(mcMoveOverlapMeter, mu, pSplitter, mu, density, numAtoms/density);
            avgP2.setDataSourceMuRoot(dsmr);
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
            dsmr1 = new DataSourceMuRoot1(mcMoveOverlapMeter, mu, pSplitter, Double.NaN, density, numAtoms/density);
            dsmr1.setMinMu(10);
            dsmr1History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            dsmr1History.setTimeDataSource(timeDataSource);
            DataPumpListener dsmr1Pump = new DataPumpListener(dsmr1, dsmr1History, numAtoms);
//                sim.integrator.getEventManager().addListener(dsmr1Pump);
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
//                sim.integrator.getEventManager().addListener(dsmrvc1Pump);
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

            getContentPane().add(simGraphic.getPanel());

            return;
    		
    	}
    }

    public static class LJParams extends ParameterBase {
        public int numAtoms = 500;
        public double density = 1.0;
        public double temperature = 0.8;
        public long steps = 1000000;
        public boolean graphics = false;
        public double mu = Double.NaN;
        public int numV = 5;
        public boolean doReweight = true;
        public double daDef = Double.NaN;
        public boolean shifted = true;
        public double rc = 3;
        public boolean fixedDaDef = false;
        public boolean ss = false;
    }
}
