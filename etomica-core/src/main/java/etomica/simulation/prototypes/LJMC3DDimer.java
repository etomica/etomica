/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeOriented;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveAtomRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P22CLJmuQAtomic;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Bar;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class LJMC3DDimer extends Simulation {

    public final PotentialComputePairGeneral potentialMaster;
    public final IntegratorMC integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public LJMC3DDimer() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public LJMC3DDimer(SimParams params) {
        super(Space3D.getInstance());

        SpeciesGeneral species = new SpeciesBuilder(space)
                .addAtom(new AtomTypeOriented(new ElementSimple("A"), Vector.of(1, 0, 0)), Vector.of(0, 0, 0))
                .build();
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density / 2);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        double Q = Debye.UNIT.toSim(-2);
        double mu = Debye.UNIT.toSim(.11);

        NeighborManagerSimple neighborManager = new NeighborManagerSimple(box);
        potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);
        P22CLJmuQAtomic p2 = new P22CLJmuQAtomic(space, 3.0058, 3.56379, 51.8037, 31.5550, mu, Q, 0.4847, 0.6461);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2);

        integrator = new IntegratorMC(potentialMaster, this.getRandom(), 1.0, box);
        integrator.setTemperature(Kelvin.UNIT.toSim(params.temperatureK));

        MCMoveAtom mcMoveMolecule = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveMolecule);

        MCMoveAtomRotate mcMoveRotate = new MCMoveAtomRotate(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveRotate);

        BoxImposePbc bip = new BoxImposePbc(space);
        bip.setBox(box);
        integrator.getEventManager().addListener(new IntegratorListenerAction(bip, params.numAtoms));

    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
        }

        LJMC3DDimer sim = new LJMC3DDimer(params);
        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperatureK);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished");

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
        sim.integrator.getEventManager().addListener(pump);

        MeterPressure meterPM = new MeterPressure(sim.box(), sim.potentialMaster);
        meterPM.setTemperature(Kelvin.UNIT.toSim(params.temperatureK));
        AccumulatorAverageFixed accp = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpP = new DataPumpListener(meterPM, accp, params.numAtoms);
        sim.integrator.getEventManager().addListener(pumpP);

        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        long t2 = System.currentTimeMillis();

        DataGroup dataPE = (DataGroup) acc.getData();
        int numAtoms = sim.getBox(0).getLeafList().size();
        double avg = dataPE.getValue(acc.AVERAGE.index) / numAtoms;
        double err = dataPE.getValue(acc.ERROR.index) / numAtoms;
        double cor = dataPE.getValue(acc.BLOCK_CORRELATION.index);

        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);

        DataGroup dataP = (DataGroup) accp.getData();
        double avgP = dataP.getValue(accp.AVERAGE.index);
        double errP = dataP.getValue(accp.ERROR.index);
        double corP = dataP.getValue(accp.BLOCK_CORRELATION.index);

        System.out.println("pressure avg: " + avgP + "  err: " + errP + "  cor: " + corP);

        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();
            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            } else {
                params.density = 0.001;
                // modify parameters here for interactive testing
            }

            LJMC3DDimer sim = new LJMC3DDimer(params);
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ((DiameterHashByType) graphic.getDisplayBox(sim.box()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 3.0058);
            ((DiameterHashByType) graphic.getDisplayBox(sim.box()).getDiameterHash()).setDiameter(sim.species().getAtomType(1), 3.56379);

            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.integrator.getEventManager().addListener(pumpPE);

            MeterPressure meterPM = new MeterPressure(sim.box(), sim.potentialMaster);
            AccumulatorAverageCollapsing accp = new AccumulatorAverageCollapsing();
            AccumulatorHistory accph = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataFork df = new DataFork();
            df.addDataSink(accp);
            df.addDataSink(accph);
            DataPumpListener pumpP = new DataPumpListener(meterPM, df, params.numAtoms);
            sim.integrator.getEventManager().addListener(pumpP);

            DisplayPlot historyPressure = new DisplayPlot();
            accph.addDataSink(historyPressure.getDataSet().makeDataSink());
            historyPressure.setLabel("Pressure");
            historyPressure.getPlot().setYLabel("Pressure (bar)");
            historyPressure.setUnit(Bar.UNIT);
            graphic.add(historyPressure);

            DisplayPlotXChart historyPressure2 = new DisplayPlotXChart();
            accph.addDataSink(historyPressure2.makeSink("Pressure"));
            historyPressure2.setLabel("Pressure");
            historyPressure2.setYLabel("Pressure (bar)");
            historyPressure2.setDefaultUnit(Bar.UNIT);
            graphic.add(historyPressure2);

            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");
            graphic.add(historyPE);

            // TODO
            DisplayPlotXChart historyPE2 = new DisplayPlotXChart();
            accPE.addDataSink(historyPE2.makeSink("PE"));
            historyPE2.setLabel("PE (XChart)");
            historyPE2.getSeries("PE").setLabel("PE");
            graphic.add(historyPE2);


            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(accp);
            displayP.setUnit(Bar.UNIT);
            displayP.setLabel("Pressure (bar)");
            graphic.add(displayP);

            graphic.makeAndDisplayFrame();

        }
    }

    public static class SimParams extends ParameterBase {
        public long steps = 1000000;
        public double density = 0.001;
        public double temperatureK = 300;
        public int numAtoms = 50;
    }
}
