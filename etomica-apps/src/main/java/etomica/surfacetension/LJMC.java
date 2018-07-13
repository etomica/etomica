/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.surfacetension;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureTensor;
import etomica.data.meter.MeterProfileByVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtomInRegion;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.List;

public class LJMC extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final SpeciesSpheresMono species;
    public final Box box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtomInRegion mcMoveAtom, mcMoveAtomBigStep;

    public LJMC(Space _space, int numAtoms, double temperature, double aspectRatio) {
        super(_space);

        //species
        species = new SpeciesSpheresMono(this, space);//index 1
        species.setIsDynamic(true);
        addSpecies(species);

        double rc = 4.0;
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);
        potentialMaster.lrcMaster().setEnabled(false);

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, random, 1.0, box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //instantiate several potentials for selection in combo-box
        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
        Vector dim = space.makeVector();
        box.getBoundary().setBoxSize(dim);

        mcMoveAtom = new MCMoveAtomInRegion(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtomBigStep = new MCMoveAtomInRegion(random, potentialMaster, space);
//        mcMoveAtomBigStep.setStepSize(6);
//        mcMoveAtomBigStep.setStepSizeMin(5);
        integrator.getMoveManager().addMCMove(mcMoveAtomBigStep);
        integrator.getMoveManager().setFrequency(mcMoveAtomBigStep, 0.2);

        //new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        double Tc = 1.31;
        double rhoc = 0.314;
        double T = temperature;
        double rhoV = rhoc - 0.477 * Math.pow(Tc - T, 1.0 / 3.0) + 0.05333 * (Tc - T) + 0.1261 * Math.pow(Tc - T, 1.5);
        double rhoL = rhoc + 0.477 * Math.pow(Tc - T, 1.0 / 3.0) + 0.2124 * (Tc - T) - 0.01151 * Math.pow(Tc - T, 1.5);
        double xL = 1.0 / (1.0 + Math.pow(rhoL / rhoV, 1.0 / 3.0));
        // V = NV/rhoV + NL/rhoL
        // rho = 0.75 rhoV + 0.25 rhoL
        double rho = (1 - xL) * rhoV + xL * rhoL;
        double V = numAtoms / rho;

        // Lx * Lyz * Lyz = V
        // Lx = Lyz * 4
        // 4 * Lyz^3 = V
        // Lyz = (V/4)^(1/3)

        double Lyz = Math.pow(V / aspectRatio, 1.0 / 3.0);
        if (Lyz < 2 * rc) {
            throw new RuntimeException("cutoff(" + rc + ") too large for box size (" + Lyz + ")");
        }
        double Lx = Lyz * aspectRatio;

        mcMoveAtom.setXRange(-0.6 * Lx * xL, 0.6 * Lx * xL, 5);
        mcMoveAtomBigStep.setXRange(0.4 * Lx * xL, -0.4 * Lx * xL, 50);
        mcMoveAtomBigStep.setStepSize(6);
        mcMoveAtomBigStep.setStepSizeMin(5);
        ((MCMoveStepTracker) mcMoveAtomBigStep.getTracker()).setAcceptanceTarget(0.2);
//        ((MCMoveStepTracker)mcMoveAtom.getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)mcMoveAtomBigStep.getTracker()).setNoisyAdjustment(true);


        //initialize box to just liquid size, put atoms on FCC lattice
        int nL = (int) ((Lx * xL * Lyz * Lyz) * rhoL);
        box.setNMolecules(species, nL);
        box.getBoundary().setBoxSize(Vector.of(new double[]{Lx * xL, Lyz, Lyz}));
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        // then expand box size

        Box boxV = new Box(space);
        addBox(boxV);
        boxV.setNMolecules(species, numAtoms - nL);
        boxV.getBoundary().setBoxSize(Vector.of(new double[]{Lx * (1 - xL), Lyz, Lyz}));
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(boxV);

        box.setNMolecules(species, numAtoms);
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < numAtoms - nL; i++) {
            Vector vx = boxV.getLeafList().get(i).getPosition();
            Vector x = atoms.get(nL + i).getPosition();
            x.E(vx);
            if (vx.getX(0) > 0) {
                x.setX(0, vx.getX(0) + 0.5 * Lx * xL);
            } else {
                x.setX(0, vx.getX(0) - 0.5 * Lx * xL);
            }
        }

        box.getBoundary().setBoxSize(Vector.of(new double[]{Lx, Lyz, Lyz}));

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) {
        LJMCParams params = new LJMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
        }
        Space space = Space3D.getInstance();
        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double aspectRatio = params.aspectRatio;
        

        final LJMC sim = new LJMC(space, numAtoms, temperature, aspectRatio);
        
        MeterPressureTensor meterPTensor = new MeterPressureTensor(sim.potentialMaster, space);
        meterPTensor.setBox(sim.box);
        meterPTensor.setTemperature(temperature);
        
        if (true) {
            SimulationGraphic ljmcGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            SimulationGraphic.makeAndDisplayFrame(ljmcGraphic.getPanel(), "LJ surface tension");
            
            IAction recenterAction = new ActionRecenter(sim);
            IntegratorListenerAction recenterActionListener = new IntegratorListenerAction(recenterAction);
            sim.integrator.getEventManager().addListener(recenterActionListener);
            recenterActionListener.setInterval(100);
            
            
            List<DataPump> dataStreamPumps = ljmcGraphic.getController().getDataStreamPumps();
            
            MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(space);
            densityProfileMeter.setBox(sim.box);
            densityProfileMeter.getXDataSource().setNValues(400);
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            meterNMolecules.setSpecies(sim.species);
            densityProfileMeter.setDataSource(meterNMolecules);
            DataFork profileFork = new DataFork();
            AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
            densityProfileAvg.setPushInterval(10);
            DataPump profilePump = new DataPump(densityProfileMeter, profileFork);
            profileFork.addDataSink(densityProfileAvg);
            DataDump profileDump = new DataDump();
            densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            IntegratorListenerAction profilePumpListener = new IntegratorListenerAction(profilePump);
            sim.integrator.getEventManager().addListener(profilePumpListener);
            profilePumpListener.setInterval(numAtoms);
            dataStreamPumps.add(profilePump);

            final FitTanh fitTanh = new FitTanh();
            densityProfileAvg.addDataSink(fitTanh, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});

            DisplayPlot profilePlot = new DisplayPlot();
            densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            fitTanh.setDataSink(profilePlot.getDataSet().makeDataSink());
            profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "density");
            profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag(), fitTanh.getTag()}, "fit");
            profilePlot.setDoLegend(true);
            profilePlot.setLabel("Density");
            ljmcGraphic.add(profilePlot);
            
            
            final FitTanh fitTanh0 = new FitTanh();
            profileFork.addDataSink(fitTanh0);

            DisplayPlot profile0Plot = new DisplayPlot();
            profileFork.addDataSink(profile0Plot.getDataSet().makeDataSink());
            fitTanh0.setDataSink(profile0Plot.getDataSet().makeDataSink());
            profile0Plot.setLegend(new DataTag[]{densityProfileMeter.getTag()}, "density");
            profile0Plot.setLegend(new DataTag[]{fitTanh0.getTag()}, "fit");
            profile0Plot.setDoLegend(true);
            profile0Plot.setLabel("Density0");
            ljmcGraphic.add(profile0Plot);
                        
            
            
            
            DataProcessor surfaceTension = new DataProcessorSurfaceTension();
            
            DataPumpListener tensionPump = new DataPumpListener(meterPTensor, surfaceTension, numAtoms);
            DataFork tensionFork = new DataFork();
            surfaceTension.setDataSink(tensionFork);
            AccumulatorAverageCollapsing surfaceTensionAvg = new AccumulatorAverageCollapsing();
            tensionFork.addDataSink(surfaceTensionAvg);
            DisplayTextBoxesCAE tensionBox = new DisplayTextBoxesCAE();
            tensionBox.setAccumulator(surfaceTensionAvg);
            tensionBox.setLabel("Surface Tension");
            tensionBox.setPrecision(4);
            dataStreamPumps.add(tensionPump);
            sim.integrator.getEventManager().addListener(tensionPump);
            ljmcGraphic.add(tensionBox);
            
            SurfaceTensionMapped stMap = new SurfaceTensionMapped(space, sim.box, sim.species, sim.potentialMaster);
            tensionFork.addDataSink(stMap);
            DataFork mappedFork = new DataFork();
            stMap.setDataSink(mappedFork);
            AccumulatorAverageCollapsing surfaceTensionMappedAvg = new AccumulatorAverageCollapsing();
            mappedFork.addDataSink(surfaceTensionMappedAvg);
            DisplayTextBoxesCAE mappedTensionBox = new DisplayTextBoxesCAE();
            mappedTensionBox.setAccumulator(surfaceTensionMappedAvg);
            mappedTensionBox.setLabel("Mapped Surface Tension");
            mappedTensionBox.setPrecision(4);
            ljmcGraphic.add(mappedTensionBox);
            
            AccumulatorHistory stHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
            stHistory.setTimeDataSource(stepCounter);
            tensionFork.addDataSink(stHistory);
            AccumulatorHistory stMappedHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            stMappedHistory.setTimeDataSource(stepCounter);
            mappedFork.addDataSink(stMappedHistory);
            
            DisplayPlot stPlot = new DisplayPlot();
            stHistory.setDataSink(stPlot.getDataSet().makeDataSink());
            stPlot.setLegend(new DataTag[]{stHistory.getTag()}, "conv");
            stMappedHistory.setDataSink(stPlot.getDataSet().makeDataSink());
            stPlot.setLegend(new DataTag[]{stMappedHistory.getTag()}, "mapped");
            stPlot.setLabel("History");
            ljmcGraphic.add(stPlot);
            
            return;
        }
        
        sim.getController().actionPerformed();

    }
    
    public static class LJMCParams extends ParameterBase {
        public int numAtoms = 1000;
        public double temperature = 1.1;
        public double aspectRatio = 6;
    }
}
