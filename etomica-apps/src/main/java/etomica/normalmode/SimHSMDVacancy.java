/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSplitter.IDataSinkFactory;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHard;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.integrator.mcmove.MCMoveIDBiasAction;
import etomica.integrator.mcmove.MCMoveInsertDeleteLatticeVacancy;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.DataSourceMuRoot.DataSourceMuRootVacancyConcentration;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.statmech.HardSphereSolid;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimHSMDVacancy extends Simulation {
    
    public final PotentialMasterList potentialMasterList;
    public final ActivityIntegrate ai;
    public IntegratorHardMDMC integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2HardSphere potential;
    public IntegratorMC integratorMC;
    public MCMoveVolume mcMoveVolume;
    public MCMoveInsertDeleteLatticeVacancy mcMoveID;
    

    public SimHSMDVacancy(final int numAtoms, double density, double tStep, int hybridInterval, final int numV, final double mu) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numAtoms);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        PrimitiveCubic primitive = new PrimitiveCubic(space, n * L);
        int[] nCells = new int[]{n, n, n};
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        BasisBigCell basis = new BasisBigCell(space, basisFCC, nCells);

        box.setBoundary(boundary);


        double nbrRange = 1.7;
        potentialMasterList = new PotentialMasterList(this, nbrRange, space);
        potentialMasterList.setCellRange(2);
        double sigma = 1.0;
        integrator = new IntegratorHardMDMC(this, potentialMasterList, box);
        integrator.setTimeStep(tStep);
        integrator.setIsothermal(true);
        integrator.setTemperature(1);
        integrator.setThermostatNoDrift(true);
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        if (nbrRange > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("rcShort is too large");
        }

        double nbr1 = L / Math.sqrt(2);
        double y = 1.25 * nbr1; //nbr1+(L-nbr1)*0.6+0.06;

        potential = new P2HardSphere(space, y, false);
        AtomType leafType = species.getLeafType();

        potentialMasterList.addPotential(potential, new AtomType[]{leafType, leafType});
        integrator.setThermostat(ThermostatType.HYBRID_MC);
        integrator.setThermostatInterval(hybridInterval);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
        // we just needed this to initial coordinates, to keep track of lattice positions of atoms
        coordinateDefinition = null;

        integrator.getEventManager().addListener(potentialMasterList.getNeighborManager(box));
        //potentialMasterList.reset();

        integratorMC = new IntegratorMC(potentialMasterList, random, 1, box);
//        mcMoveID = new MCMoveInsertDelete1(potentialMasterList, random, space, fixedN, maxDN);
//        mcMoveID = new MCMoveInsertDelete2(potentialMasterList, random, space, integrator, potentialTruncatedForceShifted, maxDN);
//        mcMoveID = new MCMoveInsertDelete3(potentialMasterList, random, space, integrator, L*Math.sqrt(0.9), fixedN, maxDN);
//        mcMoveID = new MCMoveInsertDelete4(potentialMasterList, random, space, integrator, maxDN);
        mcMoveID = new MCMoveInsertDeleteLatticeVacancy(potentialMasterList, random, space, integrator, y, numAtoms, numV);
//        integrator.setMCMoveID(mcMoveID);
//        integratorMC.getMoveEventManager().addListener(integrator);
        double x = (nbr1 - 1) / 4.0;
        mcMoveID.setMaxInsertDistance(x);
        mcMoveID.makeFccVectors(nbr1);
        mcMoveID.setMu(mu);
        mcMoveID.setSpecies(species);
        integratorMC.getMoveManager().addMCMove(mcMoveID);
        integratorMC.getMoveEventManager().addListener(mcMoveID);
        integrator.setIntegratorMC(integratorMC, numAtoms);
        potentialMasterList.reset();

        // we set the potential range high so that the MC move could find its neighbors.
        // now set potential range back to its appropriate value.
        potential.setCollisionDiameter(sigma);
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

        final SimHSMDVacancy sim = new SimHSMDVacancy(numAtoms, density, tStep, hybridInterval, numV, mu);


        final int biasInterval = 1000;

        final MCMoveOverlapListener mcMoveOverlapMeter = new MCMoveOverlapListener(sim.mcMoveID, 11, daDef, numAtoms, 2);
        mcMoveOverlapMeter.setTemperature(1);
        sim.integratorMC.getMoveEventManager().addListener(mcMoveOverlapMeter);

        final MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
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

        final MCMoveIDBiasAction mcMoveBiasAction;
        final IntegratorListenerAction mcMoveBiasListener;

        if (params.doReweight) {
            mcMoveBiasAction = new MCMoveIDBiasAction(sim.integratorMC, sim.mcMoveID, numV,
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
                double rc = sim.mcMoveID.getMaxDistance();
                double rc2 = rc*rc;
                int nmax = 12;
                public Color getAtomColor(IAtom a) {
                    if (!sim.integrator.getEventManager().firingEvent() && !sim.ai.isPaused()) return new Color(1.0f, 1.0f, 1.0f);

                    Vector pi = a.getPosition();
                    NeighborListManager nbrManager = sim.potentialMasterList.getNeighborManager(sim.box);
                    Boundary boundary = sim.box.getBoundary();
                    IAtomList nbrs = null;
                    try {
                        IAtomList[] nbrsArray = nbrManager.getUpList(a);
                        if (nbrsArray==null || nbrsArray.length==0) return Color.RED;
                        nbrs = nbrsArray[0];
                    }
                    catch (NullPointerException e) {
                        // during addition
                        return Color.RED;
                    }
                    if (nbrs==null) return Color.RED;
                    int n = 0;
                    for (int j=0; j<nbrs.getAtomCount(); j++) {
                        IAtom nbrj = nbrs.getAtom(j);
                        if (nbrj == null) break;
                        dr.Ev1Mv2(pi, nbrj.getPosition());
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();
                        if (r2 < rc2) {
                            n++;
                        }
                    }
                    nbrs = nbrManager.getDownList(a)[0];
                    for (int j=0; j<nbrs.getAtomCount(); j++) {
                        IAtom nbrj = nbrs.getAtom(j);
                        if (nbrj == null) break;
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
                Vector dr = sim.space.makeVector();
                double rc = sim.mcMoveID.getMaxDistance();
                double rc2 = rc*rc;
                int nmax = 12;
                public double getDiameter(IAtom a) {
                    Vector pi = a.getPosition();
                    NeighborListManager nbrManager = sim.potentialMasterList.getNeighborManager(sim.box);
                    Boundary boundary = sim.box.getBoundary();
                    IAtomList nbrs = null;
                    try {
                        IAtomList[] nbrsArray = nbrManager.getUpList(a);
                        if (nbrsArray==null || nbrsArray.length==0) return 0;
                        nbrs = nbrsArray[0];
                    }
                    catch (NullPointerException e) {
                        // during addition
                        return 0;
                    }
                    if (nbrs==null) return 0;
                    int n = 0;
                    for (int j=0; j<nbrs.getAtomCount(); j++) {
                        IAtom nbrj = nbrs.getAtom(j);
                        if (nbrj == null) break;
                        dr.Ev1Mv2(pi, nbrj.getPosition());
                        boundary.nearestImage(dr);
                        double r2 = dr.squared();
                        if (r2 < rc2) {
                            n++;
                        }
                    }
                    nbrs = nbrManager.getDownList(a)[0];
                    for (int j=0; j<nbrs.getAtomCount(); j++) {
                        IAtom nbrj = nbrs.getAtom(j);
                        if (nbrj == null) break;
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
            simGraphic.getDisplayBox(sim.box).setDiameterHash(dh);

            
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
    
            final AccumulatorHistory nHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSourceCountTime timeDataSource = new DataSourceCountTime(sim.integrator);
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
            DataPumpListener feHistPump = new DataPumpListener(feHistogram, nPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(feHistPump);
            nPlot.setLegend(new DataTag[]{feHistogram.getTag()}, "measured");
            
            final DataSourceImposedFEHistogram feHistogramImposed = new DataSourceImposedFEHistogram(mcMoveOverlapMeter, sim.mcMoveID, mu);
            DataPumpListener feHistImposedPump = new DataPumpListener(feHistogramImposed, nPlot.getDataSet().makeDataSink(), hybridInterval);
            nPlot.setLegend(new DataTag[]{feHistogramImposed.getTag()}, "imposed");
            sim.integrator.getEventManager().addListener(feHistImposedPump);

            if (!fluid) {
                final DataSourceImposedFEFit feHistogramFit = new DataSourceImposedFEFit(mcMoveOverlapMeter, sim.mcMoveID, mu);
                DataPumpListener feHistFitPump = new DataPumpListener(feHistogramFit, nPlot.getDataSet().makeDataSink(), hybridInterval);
                nPlot.setLegend(new DataTag[]{feHistogramFit.getTag()}, "fit");
                sim.integrator.getEventManager().addListener(feHistFitPump);
            }
            
            DisplayPlot fePlot = new DisplayPlot();
            final DataSourceFE dsfe = new DataSourceFE(mcMoveOverlapMeter);
            DataPumpListener fePump = new DataPumpListener(dsfe, fePlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(fePump);
            fePlot.setLegend(new DataTag[]{dsfe.getTag()}, "FE");

            if (!fluid) {
                final DataSourceFE dsfe2 = new DataSourceFE(mcMoveOverlapMeter);
                dsfe2.setSubtractComb(true);
                DataPumpListener fe2Pump = new DataPumpListener(dsfe2, fePlot.getDataSet().makeDataSink(), hybridInterval);
                sim.integrator.getEventManager().addListener(fe2Pump);
                fePlot.setLegend(new DataTag[]{dsfe2.getTag()}, "FEdef");
                
                final DataSourceFE dsfe3 = new DataSourceFE(mcMoveOverlapMeter);
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
            
            final DataSourceAvgPressure avgP = new DataSourceAvgPressure(pSplitter, mcMoveOverlapMeter, mu);
            DataPumpListener avgPPump = new DataPumpListener(avgP, pDisplayBox);
            sim.integrator.getEventManager().addListener(avgPPump);

            final DataSourceAvgPressure2 avgP2 = new DataSourceAvgPressure2(pSplitter, mcMoveOverlapMeter, density, numAtoms/density);
            DataPumpListener avgPPump2 = new DataPumpListener(avgP2, pDisplayBox2);
            sim.integrator.getEventManager().addListener(avgPPump2);
            
            simGraphic.add(pDisplayBox);
            simGraphic.add(pDisplayBox2);

            DisplayPlot pmuPlot = new DisplayPlot();
            pmuPlot.setDoLegend(false);
            pmuPlot.setLabel("dP vs. mu");
            final DataSourcePMu pmu = new DataSourcePMu(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, false);
            DataPumpListener pmuPump = new DataPumpListener(pmu, pmuPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(pmuPump);
            
            simGraphic.add(pmuPlot);

            DisplayPlot muMuPlot = new DisplayPlot();
            muMuPlot.setDoLegend(false);
            muMuPlot.setLabel("mu vs. mu");
            final DataSourcePMu muMu = new DataSourcePMu(mcMoveOverlapMeter, 0.2, 21, mu, pSplitter, Alat, density, numAtoms/density, true);
            DataPumpListener muMuPump = new DataPumpListener(muMu, muMuPlot.getDataSet().makeDataSink(), hybridInterval);
            sim.integrator.getEventManager().addListener(muMuPump);
            
            simGraphic.add(muMuPlot);

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
                    pmu.setMu(newValue);
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
            DataPumpListener dsmrPump = new DataPumpListener(dsmr, dsmrHistory, hybridInterval);
            sim.integrator.getEventManager().addListener(dsmrPump);
            DisplayPlot dsmrPlot = new DisplayPlot();
            dsmrPlot.setLabel("mu*");
            dsmrPlot.setLegend(new DataTag[]{dsmr.getTag()},"avg");
            dsmrHistory.setDataSink(dsmrPlot.getDataSet().makeDataSink());
            simGraphic.add(dsmrPlot);

            DataSourceMuRoot1 dsmr1 = null;
            final AccumulatorHistory dsmr1History;
            if (!fluid) {
                dsmr1 = new DataSourceMuRoot1(mcMoveOverlapMeter, mu, pSplitter, Alat, density, numAtoms/density);
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
                DataSourceMuRoot1.DataSourceMuRootVacancyConcentration dsmrvc1 = dsmr1.new DataSourceMuRootVacancyConcentration();
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

        sim.ai.setMaxSteps(steps/10);
        sim.ai.actionPerformed();

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
