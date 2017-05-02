/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.*;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.chem.elements.Iron;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterStructureFactor;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.meam.P2EAM;
import etomica.meam.PotentialCalculationEnergySumEAM;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;
import etomica.util.HistoryCollapsingAverage;
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
    public IBox box;
    public P2EAM potential;
    public P1ImageHarmonic p1ImageHarmonic;
    public MCMoveAtomSwap mcMoveSwap;

    public SimFe(Crystal crystal, int numAtoms, double temperature, double density, double w, int offsetDim, int numInnerSteps) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(space, Iron.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        IVectorMutable l = space.makeVector();
        l.E(10);
        for (int i=0; i<=offsetDim; i++) {
            l.setX(i,20);
        }
        box.getBoundary().setBoxSize(l);
        
        if (crystal != HCP) {
            BoxInflate inflater = new BoxInflate(box, space);
            inflater.setTargetDensity(density);
            inflater.actionPerformed();
        }
    
        double n = 8.7932;
        double m = 8.14475;
        double eps = ElectronVolt.UNIT.toSim(0.0220225);
        double a = 3.48501;
        double C = 28.8474;
        double rc = 6;
        potential = new P2EAM(space, n, m, eps, a, C, rc, rc);
        if (crystal == HCP) {
            int nc = (int)Math.round(Math.pow(numAtoms/8, 1.0/3.0));
            if (8*nc*nc*nc != numAtoms) {
                throw new RuntimeException("Not compatible with HCP");
            }
            // V = nc
            // v = 2/density = (2^.5/rho)sqrt(8/3)*f
            //               = 4f/(rho sqrt(3))
            // f = sqrt(3)/2
            // v = sqrt(3)/2 a^3 coa
            // a = (2 v / (sqrt(3) coa))^(1/3)
            //   = (4 / (sqrt(3) rho coa))^(1/3)
            double coa = Math.sqrt(8.0/3.0);
            double ac = Math.pow(4/(Math.sqrt(3)*density*coa), 1.0/3.0);
            double cc = coa*ac;
            IVector[] boxDim = new IVector[3];
            boxDim[0] = space.makeVector(new double[]{2*nc*ac, 0, 0});
            boxDim[1] = space.makeVector(new double[]{-2*nc*ac*Math.cos(Degree.UNIT.toSim(60)), 2*nc*ac*Math.sin(Degree.UNIT.toSim(60)), 0});
            boxDim[2] = space.makeVector(new double[]{0, 0, nc*cc});
            System.out.println("a: "+ac+" c: "+cc+" nc: "+nc);
            Primitive primitive = new PrimitiveHexagonal(space, ac, cc);
            int[] nCells = new int[]{2*nc,2*nc,nc};
            BoundaryDeformableLattice boundary = new BoundaryDeformableLattice(primitive, nCells);
            boundary.setTruncationRadius(rc);
            System.out.println(Arrays.toString(nCells));
            IVector edge0 = boundary.getEdgeVector(0);
            System.out.println(Math.sqrt(edge0.squared())+" "+edge0);
            IVector edge1 = boundary.getEdgeVector(1);
            System.out.println(Math.sqrt(edge0.squared())+" "+edge1);
            IVector edge2 = boundary.getEdgeVector(2);
            System.out.println(Math.sqrt(edge2.squared())+" "+edge2);
    
            box.setBoundary(boundary);
    
            ConfigurationLattice config = new ConfigurationLattice(new LatticeHcp(space, ac),space);
            config.setBoundaryPadding(0);
            config.initializeCoordinates(box);
        }
    
        if (crystal == HCP) {
            BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, space);
            BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class);
            potentialMaster = new PotentialMasterList(this, 1.2 * rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, space), space);
        }
        else {
            potentialMaster = new PotentialMasterList(this, 1.2*rc, space);
        }
        potentialMaster.setCellRange(2);
        double sigma = 1.0;
        if (numInnerSteps==0) {
//            integrator = new IntegratorImageHarmonicMD0(potentialMaster, random, 0.001, temperature, space);
        }
        else {
            integrator = new IntegratorImageHarmonicMD(potentialMaster, random, 0.001, temperature, space);
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

        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType,leafType});

        IVectorMutable offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim)*0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w);
        potentialMaster.addPotential(p1ImageHarmonic, new IAtomType[]{leafType});

        integrator.setBox(box);

        if (numInnerSteps>0) {
            ((IntegratorImageHarmonicMD)integrator).setP1Harmonic(p1ImageHarmonic);
            ((IntegratorImageHarmonicMD)integrator).setNumInnerSteps(numInnerSteps);
        }
        else {
//            ((IntegratorImageHarmonicMD0) integrator).setP1Harmonic(p1ImageHarmonic);
        }
        p1ImageHarmonic.setZeroForce(false);

        ConfigurationLattice config = null;
        if (crystal == Crystal.FCC) {
            config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        }
        else if (crystal == Crystal.BCC) {
            config = new ConfigurationLattice(new LatticeCubicBcc(space), space);
        }
        else if (crystal == HCP) {
            // do nothing -- already done
        }
        else {
            throw new RuntimeException("Don't know how to do "+crystal);
        }
        if (config != null) config.initializeCoordinates(box);
    
        p1ImageHarmonic.findNOffset(box);
    
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        potentialMaster.getNeighborManager(box).reset();
    
        IVector boxLength = box.getBoundary().getBoxSize();
        double lMin = boxLength.getX(0);
        if (boxLength.getX(1) < lMin) lMin = boxLength.getX(1);
        if (boxLength.getX(2) < lMin) lMin = boxLength.getX(2);
        double ww = w / lMin;
        double swapDistance = 5*Math.sqrt(1.5*temperature/ww);
        if (swapDistance > lMin/4) swapDistance = lMin/4;
        if (swapDistance < 2) swapDistance = 2;
        mcMoveSwap = new MCMoveAtomSwap(random, potentialMaster, space, p1ImageHarmonic);
        mcMoveSwap.setNbrDistance(swapDistance);
        IntegratorMC integratorMC = new IntegratorMC(potentialMaster, random, temperature);
        integrator.setIntegratorMC(integratorMC, 100);
        integrator.getIntegratorMC().getMoveManager().addMCMove(mcMoveSwap);
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 1024;
            params.steps = 10000;
            params.density = 0.15;
            params.T = 6000;
            params.w = 1000;
            params.crystal = Crystal.BCC;
            params.offsetDim = 2;
            params.numInnerSteps = 100;
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
        int numInnerSteps = params.numInnerSteps;

        if (!graphics) {
            System.out.println("Running Iron MC with N="+numAtoms+" at rho="+density+" T="+temperatureK);
            System.out.println(steps+" steps");
            System.out.println("w: "+w);
            System.out.println(numInnerSteps+" inner steps");
        }

        double L = Math.pow(numAtoms/density, 1.0/3.0);
        final SimFe sim = new SimFe(crystal, numAtoms, temperature, density, w, offsetDim, numInnerSteps);

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

            AccumulatorHistory energyHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSplitter splitter = new DataSplitter();
            DataPumpListener energyPump = new DataPumpListener(dsEnergies, splitter, 1);
            sim.integrator.getEventManager().addListener(energyPump);
            splitter.setDataSink(1, energyHist);
            DisplayPlot energyPlot = new DisplayPlot();
            energyPlot.setLabel("Energy");
            energyPlot.setUnit(new CompoundUnit(new Unit[]{ElectronVolt.UNIT, new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{1,-1}));
//            energyPlot.setUnit(ElectronVolt.UNIT);
            energyHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            simGraphic.add(energyPlot);
    
            MeterKineticEnergy meterKE = new MeterKineticEnergy();
            meterKE.setBox(sim.box);
            AccumulatorHistory keHist = new AccumulatorHistory();
            DataPumpListener kePump = new DataPumpListener(meterKE, keHist, 1);
            sim.integrator.getEventManager().addListener(kePump);
//            keHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            energyPlot.setLegend(new DataTag[]{energyHist.getTag()}, "u");
//            energyPlot.setLegend(new DataTag[]{keHist.getTag()}, "ke");

            AccumulatorAverageFixed avgSfac = new AccumulatorAverageFixed(1);
            avgSfac.setPushInterval(1);
            DataPumpListener pumpSfac = new DataPumpListener(meterSfac, avgSfac, 500);
            sim.integrator.getEventManager().addListener(pumpSfac);
            DisplayPlot plotSfac = new DisplayPlot();
            avgSfac.addDataSink(plotSfac.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{avgSfac.AVERAGE});
            plotSfac.setLabel("Structure Factor");
            plotSfac.setDoDrawLines(new DataTag[]{meterSfac.getTag()},false);
            plotSfac.getPlot().setYLog(true);
            simGraphic.add(plotSfac);


            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }

        // initial conservation of energy is often poor.  use a smaller timestep
        // for a few steps to get off lattice sites
        sim.ai.setMaxSteps(steps/20);
        sim.integrator.setTimeStep(0.0002);
        sim.getController().actionPerformed();
        sim.ai.setMaxSteps(steps/10);
        sim.integrator.setTimeStep(0.001);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.ai.setMaxSteps(steps);

        System.out.println("equilibration finished ("+steps/20+"+"+steps/10+" steps)");

        long t1 = System.currentTimeMillis();

        int interval = 10;
        long blockSize = steps/100/interval;
        if (blockSize==0) blockSize = 1;

        AccumulatorAverageFixed accEnergies = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpEnergies = new DataPumpListener(dsEnergies, accEnergies, interval);
        sim.integrator.getEventManager().addListener(pumpEnergies);

        sim.getController().actionPerformed();

        System.out.println("swap acceptance: "+sim.mcMoveSwap.getTracker().acceptanceProbability());
        System.out.println("Hybrid MD/MC acceptance: "+sim.integrator.getHybridAcceptance());

        IData avgEnergies = accEnergies.getData(accEnergies.AVERAGE);
        IData errEnergies = accEnergies.getData(accEnergies.ERROR);
        IData corEnergies = accEnergies.getData(accEnergies.BLOCK_CORRELATION);

        System.out.println("spring energy: "+avgEnergies.getValue(0)/numAtoms+"   error: "+errEnergies.getValue(0)/numAtoms+"  cor: "+corEnergies.getValue(0));
        System.out.println("Fe energy: "+avgEnergies.getValue(1)/numAtoms+"   error: "+errEnergies.getValue(1)/numAtoms+"  cor: "+corEnergies.getValue(1));

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
    }

}
