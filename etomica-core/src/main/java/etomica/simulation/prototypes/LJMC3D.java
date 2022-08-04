/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class LJMC3D extends Simulation {

    public final PotentialMasterCell potentialMaster;
    public final IntegratorMC integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public LJMC3D() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public LJMC3D(SimParams params) {
        super(Space3D.getInstance());

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        double rc = 3;
        Box box = new Box(space);
        addBox(box);
        potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());

        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones();
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);

        integrator = new IntegratorMC(potentialMaster, random, params.temperature, box);

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        this.getController().addActivity(new ActivityIntegrate(integrator));
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
        }

        LJMC3D sim = new LJMC3D(params);
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
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished");

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
        sim.integrator.getEventManager().addListener(pump);

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
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();
            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            } else {
                // modify parameters here for interactive testing
            }

            LJMC3D sim = new LJMC3D(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPE.setTimeDataSource(timeSource);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.integrator.getEventManager().addListener(pumpPE);

            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");
            graphic.add(historyPE);

            graphic.makeAndDisplayFrame();

        }
    }

    public static class SimParams extends ParameterBase {
        public long steps = 1000000;
        public double density = 0.6;
        public double temperature = 1.2;
        public int numAtoms = 1000;
    }
}
