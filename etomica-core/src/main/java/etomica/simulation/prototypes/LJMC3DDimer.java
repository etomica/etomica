/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P22CLJmuQ;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresHetero;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class LJMC3DDimer extends Simulation {

    public final PotentialMaster potentialMaster;
    public final ActivityIntegrate activityIntegrate;
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
        double rc = 3;
        potentialMaster = new PotentialMaster();

        SpeciesSpheresHetero species = new SpeciesSpheresHetero(this, space, 2);
        species.setChildCount(new int[]{1, 1});
        species.setConformation(new ConformationLinear(space, 0.4847 + 0.6461));
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density / 2);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P22CLJmuQ p2 = new P22CLJmuQ(space, 3.0058, 3.56379, 51.8037, 31.5550, 0.11, -2, 0.4847 / (0.4847 + 0.6461));
        potentialMaster.addPotential(p2, new ISpecies[]{species, species});

        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(params.temperature);

        MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(potentialMaster, random, space, 1, 1);
        integrator.getMoveManager().addMCMove(mcMoveMolecule);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // modify parameters here for interactive testing
        }

        LJMC3DDimer sim = new LJMC3DDimer(params);
        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.activityIntegrate.setMaxSteps(steps / 10);
        sim.activityIntegrate.actionPerformed();
        System.out.println("equilibration finished");

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
        sim.getIntegrator().getEventManager().addListener(pump);

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getIntegrator().resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.actionPerformed();

        long t2 = System.currentTimeMillis();

        DataGroup dataPE = (DataGroup) acc.getData();
        int numAtoms = sim.getBox(0).getLeafList().size();
        double avg = dataPE.getValue(acc.AVERAGE.index) / numAtoms;
        double err = dataPE.getValue(acc.ERROR.index) / numAtoms;
        double cor = dataPE.getValue(acc.BLOCK_CORRELATION.index);

        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();
            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            }
            else {
                params.density = 0.01;
                // modify parameters here for interactive testing
            }

            LJMC3DDimer sim = new LJMC3DDimer(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ((DiameterHashByType) graphic.getDisplayBox(sim.box()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 3.0058);
            ((DiameterHashByType) graphic.getDisplayBox(sim.box()).getDiameterHash()).setDiameter(sim.species().getAtomType(1), 3.56379);

            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.getIntegrator().getEventManager().addListener(pumpPE);

            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");
            graphic.add(historyPE);

            graphic.makeAndDisplayFrame();

        }
    }

    public static class SimParams extends ParameterBase {
        public long steps = 1000000;
        public double density = 0.01;
        public double temperature = 100;
        public int numAtoms = 1000;
    }
}
