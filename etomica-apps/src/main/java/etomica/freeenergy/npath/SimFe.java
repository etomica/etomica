/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IIntegratorListener;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Iron;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterStructureFactor;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeHcp4;
import etomica.lattice.SpaceLattice;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveHCP4;
import etomica.meam.P2EAM;
import etomica.meam.PotentialCalculationEnergySumEAM;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.Arrays;

import static etomica.freeenergy.npath.SimFe.Crystal.HCP;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimFe extends Simulation {
    
    public final PotentialMasterList potentialMaster;
    public final ActivityIntegrate ai;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2EAM potential;
    public P1ImageHarmonic p1ImageHarmonic;
    public MCMoveAtomSwap mcMoveSwap;

    public SimFe(Crystal crystal, int numAtoms, double temperature, double density, double w, int offsetDim, int numInnerSteps, boolean swap) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(space, Iron.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        Primitive primitive = (crystal == HCP) ? new PrimitiveHCP4(space) : new PrimitiveCubic(space);
        Vector l = space.makeVector();
        double[] primitiveSize = primitive.getSize();
        int[] f = new int[]{10,10,10};
        if (crystal == HCP) f = new int[]{6,4,4};
        for (int i=0; i<3; i++) {
            double x = f[i]*primitiveSize[i];
            if (i<=offsetDim) x *= 2;
            l.setX(i,x);
        }
        box.getBoundary().setBoxSize(l);
        
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        double n = 8.7932;
        double m = 8.14475;
        double eps = ElectronVolt.UNIT.toSim(0.0220225);
        double a = 3.48501;
        double C = 28.8474;
        double rc = 6;
        potential = new P2EAM(space, n, m, eps, a, C, rc, rc);

        potentialMaster = new PotentialMasterList(this, 1.2*rc, space);
        potentialMaster.getNbrCellManager(box).setSuppressBoxLengthWarning(true);
        potentialMaster.setCellRange(2);

        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        Vector offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim) * 0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w, true);
        potentialMaster.addPotential(p1ImageHarmonic, new AtomType[]{leafType});

        if (numInnerSteps > 0 || swap) {
            integrator = new IntegratorImageHarmonicMD(potentialMaster, random, 0.001, temperature, space);
            ((IntegratorImageHarmonicMD) integrator).setP1Harmonic(p1ImageHarmonic);
            ((IntegratorImageHarmonicMD) integrator).setNumInnerSteps(numInnerSteps);
        } else {
            integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.001, temperature, space);
        }
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setPotentialCalculation(new PotentialCalculationEnergySumEAM(potential));
        integrator.setMeterPotentialEnergy(meterPE);
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMD.ThermostatType.HYBRID_MC);
        integrator.setThermostatInterval(10);
        integrator.setTemperature(temperature);
        integrator.getEventManager().addListener(potential.makeIntegratorListener(potentialMaster, box));
        integrator.setForceSum(new PotentialCalculationForceSum());

        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        integrator.setBox(box);

        p1ImageHarmonic.setZeroForce(false);

        SpaceLattice lat = null;
        if (crystal == Crystal.FCC) {
            lat = new LatticeCubicFcc(space);
        }
        else if (crystal == Crystal.BCC) {
            lat = new LatticeCubicBcc(space);
        }
        else if (crystal == HCP) {
            lat = new LatticeHcp4(space);
        }
        else {
            throw new RuntimeException("Don't know how to do "+crystal);
        }
        ConfigurationLattice config = new ConfigurationLattice(lat, space);
        config.initializeCoordinates(box);
    
        p1ImageHarmonic.findNOffset(box);
    
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        potentialMaster.getNeighborManager(box).reset();
    
        Vector boxLength = box.getBoundary().getBoxSize();
        double lMin = boxLength.getX(0);
        if (boxLength.getX(1) < lMin) lMin = boxLength.getX(1);
        if (boxLength.getX(2) < lMin) lMin = boxLength.getX(2);
        double ww = w / lMin;
        double swapDistance = 1.5*Math.sqrt(1.5*temperature/ww);
        if (swapDistance > lMin/4) swapDistance = lMin/4;
        if (swapDistance > rc) swapDistance = rc;
        if (swapDistance < 2) swapDistance = 2;
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(this, swapDistance, space);
        potentialMasterCell.setCellRange(2);
        potentialMasterCell.getNbrCellManager(box).assignCellAll();
        if (swap) {
            mcMoveSwap = new MCMoveAtomSwap(random, potentialMasterCell, space, p1ImageHarmonic);
            mcMoveSwap.setNbrDistance(swapDistance);
            IntegratorMC integratorMC = new IntegratorMC(potentialMaster, random, temperature);
            integrator.setIntegratorMC(integratorMC, 10 * numAtoms);
            integrator.getIntegratorMC().getMoveManager().addMCMove(mcMoveSwap);

            integrator.getIntegratorMC().getMoveEventManager().addListener(potentialMasterCell.getNbrCellManager(box).makeMCMoveListener());
            ((IntegratorImageHarmonicMD) integrator).setNeighborCellManager(potentialMasterCell.getNbrCellManager(box));
        }
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 1024;
            params.steps = 1000;
            params.density = 0.15151515151515;
            params.T = 7000;
            params.w = 0;
            params.crystal = Crystal.BCC;
            params.offsetDim = 2;
            params.numInnerSteps = 0;
            params.swap = false;
        }

        final int numAtoms = params.numAtoms;
        final double temperatureK = params.T;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final double density = params.density;
        long steps = params.steps;
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;
        Crystal crystal = params.crystal;
        boolean swap = params.swap;
        int numInnerSteps = w > 0 ? params.numInnerSteps : (swap ? 1 : 0);

        if (!graphics) {
            System.out.println("Running Iron MC with N="+numAtoms+" at rho="+density+" T="+temperatureK);
            System.out.println(steps+" steps");
            System.out.println("w: "+w);
            System.out.println(numInnerSteps+" inner steps");
        }

        double L = Math.pow(numAtoms/density, 1.0/3.0);
        final SimFe sim = new SimFe(crystal, numAtoms, temperature, density, w, offsetDim, numInnerSteps, swap);
        System.out.println(Arrays.toString(sim.getRandomSeeds()));

        DataSourceEnergies dsEnergies = new DataSourceEnergies(sim.potentialMaster);
        dsEnergies.setPotentialCalculation(new DataSourceEnergies.PotentialCalculationEnergiesEAM(sim.potential));
        dsEnergies.setBox(sim.box);
        IData u = dsEnergies.getData();
        if (!graphics) System.out.println("Fe lattice energy (eV/atom): "+ElectronVolt.UNIT.fromSim(u.getValue(1)/numAtoms));

        MeterStructureFactor meterSfac = new MeterStructureFactor(sim.space, sim.box, 8);
        if (graphics) {
            final String APP_NAME = "SimFe";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());
            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    return a.getLeafIndex() < numAtoms/2 ? Color.RED : Color.BLUE;
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            DiameterHash diameter = simGraphic.getDisplayBox(sim.box).getDiameterHash();
            ((DiameterHashByType)diameter).setDiameter(sim.species.getLeafType(),2);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            DataSourceCountTime tSource = new DataSourceCountTime(sim.integrator);
            AccumulatorHistory energyHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyHist.setTimeDataSource(tSource);
            AccumulatorHistory springHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            springHist.setTimeDataSource(tSource);
            DataSplitter splitter = new DataSplitter();
            DataPumpListener energyPump = new DataPumpListener(dsEnergies, splitter, 10);
            sim.integrator.getEventManager().addListener(energyPump);
            splitter.setDataSink(0, springHist);
            splitter.setDataSink(1, energyHist);
            DisplayPlot energyPlot = new DisplayPlot();
            energyPlot.setLabel("Fe");
            energyPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
//            energyPlot.setUnit(new CompoundUnit(new Unit[]{ElectronVolt.UNIT, new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{1,-1}));
//            energyPlot.setUnit(ElectronVolt.UNIT);
            energyHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            simGraphic.add(energyPlot);
            DisplayPlot springPlot = new DisplayPlot();
            springPlot.setLabel("spring");
            springPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
            springHist.addDataSink(springPlot.getDataSet().makeDataSink());
            simGraphic.add(springPlot);
    
            MeterKineticEnergy meterKE = new MeterKineticEnergy();
            meterKE.setBox(sim.box);
            AccumulatorHistory keHist = new AccumulatorHistory();
            DataPumpListener kePump = new DataPumpListener(meterKE, keHist, 1);
//            sim.integrator.getEventManager().addListener(kePump);
//            keHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            energyPlot.setLegend(new DataTag[]{energyHist.getTag()}, "u");
//            energyPlot.setLegend(new DataTag[]{keHist.getTag()}, "ke");

            AccumulatorHistory dudwHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            splitter.setDataSink(2, dudwHist);
            dudwHist.setTimeDataSource(tSource);
            DisplayPlot dudwPlot = new DisplayPlot();
            dudwHist.addDataSink(dudwPlot.getDataSet().makeDataSink());
            dudwPlot.setLabel("dudw");
            dudwPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
            simGraphic.add(dudwPlot);

            AccumulatorAverageFixed avgSfac = new AccumulatorAverageFixed(1);
            avgSfac.setPushInterval(1);
            DataPumpListener pumpSfac = new DataPumpListener(meterSfac, avgSfac, 500);
            sim.integrator.getEventManager().addListener(pumpSfac);
            DisplayPlot plotSfac = new DisplayPlot();
            avgSfac.addDataSink(plotSfac.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});
            plotSfac.setLabel("Structure Factor");
            plotSfac.setDoDrawLines(new DataTag[]{meterSfac.getTag()},false);
            plotSfac.getPlot().setYLog(true);
            simGraphic.add(plotSfac);
            
            sim.integrator.setTimeStep(0.0001);
            sim.integrator.getEventManager().addListener(new IIntegratorListener() {
                @Override
                public void integratorInitialized(IntegratorEvent e) {}

                @Override
                public void integratorStepStarted(IntegratorEvent e) {}

                @Override
                public void integratorStepFinished(IntegratorEvent e) {
                    if (sim.integrator.getStepCount() > 400) {
                        sim.integrator.setTimeStep(0.001);
                    }
                }
            });


            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }

        // initial conservation of energy is often poor.  use a smaller timestep
        // for a few steps to get off lattice sites
        sim.ai.setMaxSteps(steps/20);
        sim.integrator.setTimeStep(0.0001);
        sim.getController().actionPerformed();
        sim.getController().reset();
        if (sim.integrator.getHybridAcceptance() < 0.5) {
            throw new RuntimeException("hybrid acceptance "+sim.integrator.getHybridAcceptance());
        }
        sim.integrator.resetStepCount();
        sim.ai.setMaxSteps(steps/20);
        sim.integrator.setTimeStep(0.0002);
        sim.getController().actionPerformed();
        sim.getController().reset();
        if (sim.integrator.getHybridAcceptance() < 0.5) {
            throw new RuntimeException("hybrid acceptance "+sim.integrator.getHybridAcceptance());
        }
        sim.integrator.resetStepCount();
        sim.ai.setMaxSteps(steps/10);
        sim.integrator.setTimeStep(0.001);
        sim.getController().actionPerformed();
        sim.getController().reset();
        if (sim.integrator.getHybridAcceptance() < 0.5) {
            throw new RuntimeException("hybrid acceptance "+sim.integrator.getHybridAcceptance());
        }
        sim.integrator.resetStepCount();
        sim.ai.setMaxSteps(steps);

        System.out.println("equilibration finished ("+steps/20+"+"+steps/20+"+"+steps/10+" steps)");

        long t1 = System.currentTimeMillis();

        int interval = 10;
        long blockSize = steps/100/interval;
        if (blockSize==0) blockSize = 1;

        AccumulatorAverageCovariance accEnergies = new AccumulatorAverageCovariance(blockSize);
        DataFork energyFork = new DataFork();
        DataPumpListener pumpEnergies = new DataPumpListener(dsEnergies, energyFork, interval);
        energyFork.addDataSink(accEnergies);
        sim.integrator.getEventManager().addListener(pumpEnergies);
        DataLogger dataLogger = null;
        if (!swap && w == 0) {
            dataLogger = new DataLogger();
            DataSplitter splitter = new DataSplitter();
            energyFork.addDataSink(splitter);
            splitter.setDataSink(2, dataLogger);
            dataLogger.setFileName("r2hist.dat");
            DataArrayWriter writer = new DataArrayWriter();
            writer.setIncludeHeader(false);
            dataLogger.setDataSink(writer);
            dataLogger.setAppending(false);
            sim.getController().getEventManager().addListener(dataLogger);
        }

        sim.getController().actionPerformed();

        if (swap) System.out.println("swap acceptance: " + sim.mcMoveSwap.getTracker().acceptanceProbability());
        System.out.println("Hybrid MD/MC acceptance: "+sim.integrator.getHybridAcceptance());

        IData avgEnergies = accEnergies.getData(AccumulatorAverage.AVERAGE);
        IData errEnergies = accEnergies.getData(AccumulatorAverage.ERROR);
        IData corEnergies = accEnergies.getData(AccumulatorAverage.BLOCK_CORRELATION);
        IData covEnergies = accEnergies.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);

        System.out.println("spring energy: "+avgEnergies.getValue(0)/numAtoms+"   error: "+errEnergies.getValue(0)/numAtoms+"  cor: "+corEnergies.getValue(0));
        System.out.println("Fe energy: "+avgEnergies.getValue(1)/numAtoms+"   error: "+errEnergies.getValue(1)/numAtoms+"  cor: "+corEnergies.getValue(1));
        System.out.println("du/dw: "+avgEnergies.getValue(2)/numAtoms+"   error: "+errEnergies.getValue(2)/numAtoms+"  cor: "+corEnergies.getValue(2));
        double var0 = covEnergies.getValue(0*3+0);
        double var1 = covEnergies.getValue(1*3+1);
        double var2 = covEnergies.getValue(2*3+2);
        double cor01 = 0;
        if (var0*var1>0) cor01 = covEnergies.getValue(0*3+1)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(1*3+1));
        double cor02 = 0;
        if (var0*var2>0) cor02 = covEnergies.getValue(0*3+2)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(2*3+2));
        System.out.println("spring correlation: 1 "+cor01+" "+cor02);
        double cor12 = covEnergies.getValue(1*3+2)/Math.sqrt(covEnergies.getValue(1*3+1)*covEnergies.getValue(2*3+2));
        System.out.println("Fe correlation: "+cor01+" 1 "+cor12);
        System.out.println("du/dw correlation: "+cor02+" "+cor12+" 1");

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }

    enum Crystal {FCC, BCC, HCP}

    public static class LjMC3DParams extends ParameterBase {
        public int numAtoms = 500;
        public double T = 2.0;
        public double density = 0.3;
        public long steps = 100000;
        public double rc = 2.5;
        public boolean graphics = false;
        public double w = 1;
        public int offsetDim = 0;
        public int numInnerSteps = 10;
        public Crystal crystal = Crystal.FCC;
        public boolean swap = true;
    }

}
