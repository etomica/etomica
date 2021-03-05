/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomSourceRandomSpecies;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.*;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.Potential1;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomNumberGenerator;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LJMD extends Simulation {
    
    public final PotentialMasterList potentialMaster;
    public IntegratorFixedWall integrator;
    public SpeciesGeneral speciesFluid, speciesTopWall, speciesBottomWall;
    public Box box;
    public P2SoftSphericalTruncatedForceShifted pFF, pTW, pBW;
    public ConfigurationLammps config;
    public MCMoveAtomNbr mcMove;

    public LJMD(double temperature, double tStep, boolean fixedWall, double spring, double springPosition, double Psat, int hybridInterval, int mcSteps, String lammpsFile) {
        super(Space3D.getInstance());
        setRandom(new RandomNumberGenerator(2));

        speciesFluid = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("F")), true);
        addSpecies(speciesFluid);
        speciesTopWall = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("TW", fixedWall ? Double.POSITIVE_INFINITY : 1)), true);
        addSpecies(speciesTopWall);
        speciesBottomWall = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("BW", Double.POSITIVE_INFINITY)), true);
        addSpecies(speciesBottomWall);

        BoundaryRectangularSlit boundary = new BoundaryRectangularSlit(2, space);
        box = this.makeBox(boundary);

        config = new ConfigurationLammps(space, lammpsFile, speciesTopWall, speciesBottomWall, speciesFluid);
        config.setTopPadding(50);
        config.initializeCoordinates(box);

        potentialMaster = new PotentialMasterList(this, 1.2 * 5.49925, space);
        potentialMaster.setCellRange(2);
        integrator = new IntegratorFixedWall(potentialMaster, random, tStep, temperature, box);
        integrator.setIsothermal(true);
        integrator.setTemperature(temperature);
        integrator.setThermostat(hybridInterval > 0 ? ThermostatType.HYBRID_MC : ThermostatType.ANDERSEN);
        integrator.setThermostatInterval(hybridInterval > 0 ? hybridInterval : 2000);


        pFF = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.0, 1.0), 2.5);
        AtomType leafType = speciesFluid.getLeafType();
        potentialMaster.addPotential(pFF, new AtomType[]{leafType, leafType});

        pBW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.09985, 0.4), 5.49925);
        potentialMaster.addPotential(pBW, new AtomType[]{leafType, speciesBottomWall.getLeafType()});

        pTW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.5, 0.1), 1.68);
        potentialMaster.addPotential(pTW, new AtomType[]{leafType, speciesTopWall.getLeafType()});

        if (!fixedWall) {
            FixedWall fixedWallListener = new FixedWall(space, box, integrator.getAgentManager(), speciesTopWall);
            integrator.setFixedWall(fixedWallListener);
            int nWall = box.getMoleculeList(speciesTopWall).size();
            double Lxy = boundary.getBoxSize().getX(0);
            P1Wall p1Wall = new P1Wall(space, spring / nWall, springPosition + config.getShift().getX(2), Psat * Lxy * Lxy / nWall);
            potentialMaster.addPotential(p1Wall, new AtomType[]{speciesTopWall.getLeafType()});
        }

        if (mcSteps > 0 && hybridInterval > 0) {
            IntegratorMC integratorMC = new IntegratorMC(this.getRandom(), potentialMaster, box);
            integratorMC.setTemperature(temperature);
            mcMove = new MCMoveAtomNbr(random, potentialMaster, space);
            mcMove.setAtomSource(new AtomSourceRandomSpecies(getRandom(), speciesFluid));
            mcMove.setStepSize(0.5 * config.getLxy());
            ((MCMoveStepTracker) mcMove.getTracker()).setTunable(false);

            integratorMC.getMoveManager().addMCMove(mcMove);

            integrator.setIntegratorMC(integratorMC, mcSteps);

            Potential1 p1F = new Potential1(space) {

                public double energy(IAtomList atoms) {
                    double pz = atoms.get(0).getPosition().getX(2);
                    double zMin = -0.5 * boundary.getBoxSize().getX(2);
                    double zMax = box.getMoleculeList(speciesTopWall).get(0).getChildList().get(0).getPosition().getX(2);
                    return (pz < zMin || pz > zMax) ? Double.POSITIVE_INFINITY : 0;
                }
            };
            potentialMaster.addPotential(p1F, new AtomType[]{leafType});
        }

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }
    
    public static void main(String[] args) {

        LJMDParams params = new LJMDParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.lammpsFile = "eq.data";
            params.steps = 10000;
            params.T = 0.8;
            params.fixedWall = false;
            params.springPosition = 76;
            params.hybridInterval = 0;
            params.tStep = 0.0025;
            params.mcSteps = 0;
        }

        final double temperature = params.T;
        long steps = params.steps;
        double tStep = params.tStep;
        boolean fixedWall = params.fixedWall;
        double spring = params.spring;
        double springPosition = params.springPosition;
        double Psat = params.Psat;
        int hybridInterval = params.hybridInterval;
        int mcSteps = params.mcSteps;
        int foo = params.dataInterval;
        String lammpsFile = params.lammpsFile;
        boolean graphics = params.graphics;

        if (!graphics) System.out.println(steps+" steps");

        final LJMD sim = new LJMD(temperature, tStep, fixedWall, spring, springPosition, Psat, hybridInterval, mcSteps, lammpsFile);

        if (hybridInterval > 0) {
            foo = (foo/hybridInterval)*hybridInterval;
            if (foo == 0) foo = hybridInterval;
            if (foo != params.dataInterval) System.out.println("dataInterval => "+foo);
        }
        final int dataInterval = foo;
        
        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        DataFork forkPE = new DataFork();
        DataPumpListener pumpPE = new DataPumpListener(meterPE, forkPE, dataInterval);
        sim.integrator.getEventManager().addListener(pumpPE);

        final MeterWallForce meterWF = new MeterWallForce(sim.space, sim.potentialMaster, sim.box, sim.speciesTopWall);
        DataFork forkWF = new DataFork();
        DataPumpListener pumpWF = new DataPumpListener(meterWF, forkWF, dataInterval);
        sim.integrator.getEventManager().addListener(pumpWF);
        
        double zShift = sim.config.getShift().getX(2);
        final MeterWallPosition meterWP = new MeterWallPosition(sim.box, sim.speciesTopWall, zShift);
        DataFork forkWP = new DataFork();
        DataPumpListener pumpWP = new DataPumpListener(meterWP, forkWP, dataInterval);
        sim.integrator.getEventManager().addListener(pumpWP);
        
        sim.integrator.reset();
        double u = meterPE.getDataAsScalar();
        System.out.println("Initial Potential energy: "+u);
        System.out.println("Initial Wall force: "+meterWF.getDataAsScalar());
        System.out.println("Initial Wall position: "+meterWP.getDataAsScalar());
        
        MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(sim.space);
        densityProfileMeter.setProfileDim(2);
        densityProfileMeter.setBox(sim.box);
        densityProfileMeter.setSpecies(sim.speciesFluid);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, dataInterval);
        DataDump profileDump = new DataDump();
        densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        sim.integrator.getEventManager().addListener(profilePump);
        densityProfileAvg.setPushInterval(1);

        AccumulatorHistogram histogramWP = new AccumulatorHistogram(new HistogramExpanding(0.2));
        forkWP.addDataSink(histogramWP);

        if (graphics) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            final String APP_NAME = "LJMD";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

            List<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            dataStreamPumps.add(profilePump);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHashByType dh = (DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dh.setDiameter(sim.speciesFluid.getLeafType(), 1.0);
            dh.setDiameter(sim.speciesBottomWall.getLeafType(), 1.09885);
            dh.setDiameter(sim.speciesTopWall.getLeafType(), 1.5);
            
            DataSourceCountTime dsSteps = new DataSourceCountTime(sim.integrator);

            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(dsSteps);
            forkPE.addDataSink(historyPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.setDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            simGraphic.add(plotPE);
            
            AccumulatorHistory historyWF = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyWF.setTimeDataSource(dsSteps);
            forkWF.addDataSink(historyWF);
            DisplayPlot plotWF = new DisplayPlot();
            historyWF.setDataSink(plotWF.getDataSet().makeDataSink());
            plotWF.setLabel("Force");
            simGraphic.add(plotWF);

            if (!fixedWall) {
                AccumulatorHistory historyWP = new AccumulatorHistory(new HistoryCollapsingAverage());
                historyWP.setTimeDataSource(dsSteps);
                forkWP.addDataSink(historyWP);
                DisplayPlot plotWP = new DisplayPlot();
                historyWP.setDataSink(plotWP.getDataSet().makeDataSink());
                plotWP.setLabel("Position");
                simGraphic.add(plotWP);
            }
            
            DisplayPlot profilePlot = new DisplayPlot();
            densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            profilePlot.setLabel("density");
            simGraphic.add(profilePlot);

            DisplayPlot wallHistogramPlot = new DisplayPlot();
            histogramWP.addDataSink(wallHistogramPlot.getDataSet().makeDataSink());
            wallHistogramPlot.setLabel("wall");
            simGraphic.add(wallHistogramPlot);

            return;
        }

        long bs = steps/(100*dataInterval);
        if (bs==0) bs=dataInterval;
        AccumulatorAverageFixed accPE = new AccumulatorAverageFixed(bs);
        forkPE.addDataSink(accPE);
        AccumulatorAverageFixed accWF = new AccumulatorAverageFixed(bs);
        forkWF.addDataSink(accWF);
        
        final FileWriter fw;
        try {
            fw = new FileWriter("history.dat");
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        sim.integrator.getEventManager().addListener(new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                long step = sim.integrator.getStepCount();
                if (step % dataInterval != 0) return;
                writeIt();
            }
            protected void writeIt() {
                long step = sim.integrator.getStepCount();
                double u = meterPE.getDataAsScalar();
                double wf = meterWF.getDataAsScalar();
                double wp = meterWP.getDataAsScalar();
                try {
                    fw.write(step+" "+u+" "+wf+" "+wp+"\n");
                }
                catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }
            
            public void integratorInitialized(IntegratorEvent e) {
                writeIt();
            }
        });

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        try {
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        
        u = meterPE.getDataAsScalar();
        System.out.println("Final Potential energy: "+u);
        System.out.println("Final Wall force: "+meterWF.getDataAsScalar());
        System.out.println("Final Wall position: "+meterWP.getDataAsScalar());
        if (hybridInterval>0 && mcSteps>0) System.out.println("MC acceptance "+sim.mcMove.getTracker().acceptanceProbability());
        
        FileWriter fwProfile;
        try {
            fwProfile = new FileWriter("density.dat");
            IData profileAvg = densityProfileAvg.getData(densityProfileAvg.AVERAGE);
            IData xProfile = ((DataInfoFunction) ((DataInfoGroup) densityProfileAvg.getDataInfo()).getSubDataInfo(densityProfileAvg.AVERAGE.index)).getXDataSource().getIndependentData(0);
            for  (int i=0; i<xProfile.getLength(); i++) {
                fwProfile.write((xProfile.getValue(i)-sim.config.getShift().getX(2))+" "+profileAvg.getValue(i)+"\n");
            }
            fwProfile.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        
        
        if (fixedWall) {
            double avgPE = accPE.getData().getValue(accPE.AVERAGE.index);
            double errPE = accPE.getData().getValue(accPE.ERROR.index);
            double corPE = accPE.getData().getValue(accPE.BLOCK_CORRELATION.index);
            double avgWF = accWF.getData().getValue(accWF.AVERAGE.index);
            double errWF = accWF.getData().getValue(accWF.ERROR.index);
            double corWF = accWF.getData().getValue(accWF.BLOCK_CORRELATION.index);
            
            if (steps>100*dataInterval) {
                System.out.print(String.format("Average potential energy: %25.15e %10.4e % 5.3f\n",avgPE,errPE,corPE));
                System.out.print(String.format("Average wall force: %25.15e %10.4e % 5.3f\n",avgWF,errWF,corWF));
            }
            else {
                System.out.println("Average potential energy: "+avgPE);
                System.out.println("Average wall force: "+avgWF);
            }
        }
        
        WriteConfigurationInterfacial configWriter = new WriteConfigurationInterfacial(sim.space);
        configWriter.setSpecies(sim.speciesFluid);
        Vector unshift = sim.space.makeVector();
        unshift.Ea1Tv1(-1, sim.config.getShift());
        configWriter.setShift(unshift);
        configWriter.setBox(sim.box);
        configWriter.setFileName("xyz_000.dat");
        configWriter.actionPerformed();
    }
    
    public static class LJMDParams extends ParameterBase {
        public double T = 2.0;
        public long steps = 10000;
        public double tStep = 0.0025;
        public boolean graphics = false;
        public String lammpsFile = "";
        public boolean fixedWall = true;
        public double spring = 0.3;
        public double springPosition = 70;
        public double Psat = 0.030251;
        public int hybridInterval = 0;
        public int mcSteps = 0;
        public int dataInterval = 10;
    }
}
