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
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterProfileByVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.Potential1;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LJMC extends Simulation {
    
    public final PotentialMasterCell potentialMasterCell;
    public final ActivityIntegrate ai;
    public IntegratorMC integrator;
    public SpeciesSpheresMono speciesFluid, speciesTopWall, speciesBottomWall;
    public Box box;
    public P2SoftSphericalTruncatedForceShifted pFF, pTW, pBW;
    public ConfigurationLammps config;


    public LJMC(double temperature, String lammpsFile) {
        super(Space3D.getInstance());
        BoundaryRectangularSlit boundary = new BoundaryRectangularSlit(2, space);
        box = new Box(boundary, space);
        addBox(box);
        
        speciesFluid = new SpeciesSpheresMono(space, new ElementSimple("F"));
        addSpecies(speciesFluid);
        speciesTopWall = new SpeciesSpheresMono(space, new ElementSimple("TW"));
        addSpecies(speciesTopWall);
        speciesBottomWall = new SpeciesSpheresMono(space, new ElementSimple("BW"));
        addSpecies(speciesBottomWall);
        
        config = new ConfigurationLammps(space, lammpsFile, speciesTopWall, speciesBottomWall, speciesFluid);
        config.setTopPadding(10);
        config.initializeCoordinates(box);
        double Lxy = config.getLxy();
        
        potentialMasterCell = new PotentialMasterCell(this, 5.49925, space);
        potentialMasterCell.setCellRange(2);
        integrator = new IntegratorMC(this, potentialMasterCell);
        integrator.setTemperature(temperature);
        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMasterCell, space);
        mcMoveAtom.setAtomSource(new AtomSourceRandomSpecies(getRandom(), speciesFluid));
        MCMoveAtom mcMoveAtomBig = new MCMoveAtom(random, potentialMasterCell, space);
        mcMoveAtomBig.setAtomSource(new AtomSourceRandomSpecies(getRandom(), speciesFluid));
        mcMoveAtomBig.setStepSize(0.5*Lxy);
        ((MCMoveStepTracker)mcMoveAtomBig.getTracker()).setTunable(false);
        
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().addMCMove(mcMoveAtomBig);
        
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        pFF = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.0, 1.0), 2.5);
        AtomType leafType = speciesFluid.getLeafType();
        potentialMasterCell.addPotential(pFF, new AtomType[]{leafType, leafType});

        pBW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.09985, 0.4), 5.49925);
        potentialMasterCell.addPotential(pBW, new AtomType[]{leafType, speciesBottomWall.getLeafType()});
        
        pTW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.5, 0.1), 1.68);
        potentialMasterCell.addPotential(pTW, new AtomType[]{leafType, speciesTopWall.getLeafType()});

        Potential1 p1F = new Potential1(space) {
            
            public double energy(IAtomList atoms) {
                double pz = atoms.getAtom(0).getPosition().getX(2);
                double zMin = -0.5*boundary.getBoxSize().getX(2);
                double zMax = box.getMoleculeList(speciesTopWall).getMolecule(0).getChildList().getAtom(0).getPosition().getX(2);
                return (pz < zMin || pz > zMax) ? Double.POSITIVE_INFINITY : 0;
            }
        };
        potentialMasterCell.addPotential(p1F, new AtomType[]{leafType});
        
        integrator.setBox(box);

        integrator.getMoveEventManager().addListener(potentialMasterCell.getNbrCellManager(box).makeMCMoveListener());
		
        potentialMasterCell.getNbrCellManager(box).assignCellAll();
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.lammpsFile = "eq.data";
            params.steps = 10000;
            params.T = 0.8;
        }

        final double temperature = params.T;
        long steps = params.steps;
        boolean graphics = params.graphics;
        String lammpsFile = params.lammpsFile;

        if (!graphics) {
            System.out.println("Running MC with T="+temperature);
            System.out.println(steps+" steps");
        }

        final LJMC sim = new LJMC(temperature, lammpsFile);

        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator();
        meterPE.setIntegrator(sim.integrator);
        DataFork forkPE = new DataFork();
        DataPumpListener pumpPE = new DataPumpListener(meterPE, forkPE, 10);
        sim.integrator.getEventManager().addListener(pumpPE);

        MeterWallForce meterWF = new MeterWallForce(sim.space, sim.potentialMasterCell, sim.box, sim.speciesTopWall);
        DataFork forkWF = new DataFork();
        DataPumpListener pumpWF = new DataPumpListener(meterWF, forkWF, 1000);
        sim.integrator.getEventManager().addListener(pumpWF);
        
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMasterCell);
        meterPE2.setBox(sim.box);
        double u = meterPE2.getDataAsScalar();
        System.out.println("Potential energy: "+u);
        System.out.println("Wall force: "+meterWF.getDataAsScalar());
        
        MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(sim.space);
        densityProfileMeter.setProfileDim(2);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.speciesFluid);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataPump profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 1000);
        DataDump profileDump = new DataDump();
        densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});
        IntegratorListenerAction profilePumpListener = new IntegratorListenerAction(profilePump);
        sim.integrator.getEventManager().addListener(profilePumpListener);
        profilePumpListener.setInterval(10);
        
        
        if (graphics) {
            final String APP_NAME = "LJMC";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

            List<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            dataStreamPumps.add(profilePump);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHashByType dh = (DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dh.setDiameter(sim.speciesFluid.getLeafType(), 1.0);
            dh.setDiameter(sim.speciesBottomWall.getLeafType(), 1.09885);
            dh.setDiameter(sim.speciesTopWall.getLeafType(), 1.5);
            
            DataSourceCountSteps dsSteps = new DataSourceCountSteps(sim.integrator);

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
            

            DisplayPlot profilePlot = new DisplayPlot();
            densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});
            profilePlot.setLabel("density");
            simGraphic.add(profilePlot);
            
            return;
        }

        long bs = steps/100000;
        if (bs==0) bs=1;
        AccumulatorAverageFixed accPE = new AccumulatorAverageFixed(bs);
        forkPE.addDataSink(accPE);
        AccumulatorAverageFixed accWF = new AccumulatorAverageFixed(bs);
        forkWF.addDataSink(accWF);

        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        
        u = meterPE2.getDataAsScalar();
        System.out.println("Potential energy: "+u);
        System.out.println("Wall force: "+meterWF.getDataAsScalar());

        double avgPE = accPE.getData().getValue(AccumulatorAverage.AVERAGE.index);
        double errPE = accPE.getData().getValue(AccumulatorAverage.ERROR.index);
        double corPE = accPE.getData().getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        double avgWF = accWF.getData().getValue(AccumulatorAverage.AVERAGE.index);
        double errWF = accWF.getData().getValue(AccumulatorAverage.ERROR.index);
        double corWF = accWF.getData().getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        
        if (steps>100000) {
            System.out.println(String.format("Average potential energy: %25.15e %10.4e % 5.3f\n",avgPE,errPE,corPE));
            System.out.println(String.format("Average wall force: %25.15e %10.4e % 5.3f\n",avgWF,errWF,corWF));
        }
        else {
            System.out.println("Average potential energy: "+avgPE);
            System.out.println("Average wall force: "+avgWF);
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

    public static class LjMC3DParams extends ParameterBase {
        public double T = 2.0;
        public long steps = 100000;
        public boolean graphics = false;
        public String lammpsFile = "";
    }
}
