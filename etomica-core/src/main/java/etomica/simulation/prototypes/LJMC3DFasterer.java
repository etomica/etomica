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
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCellFasterer;
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

public class LJMC3DFasterer extends Simulation {

    public final PotentialMasterCellFasterer potentialMaster;
    public final IntegratorMCFasterer integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public LJMC3DFasterer() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public LJMC3DFasterer(SimParams params) {
        super(Space3D.getInstance());

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        double rc = 3;
        Box box = new Box(space);
        addBox(box);
        potentialMaster = new PotentialMasterCellFasterer(getSpeciesManager(), box, 2, BondingInfo.noBonding());

        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones(space);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);

        integrator = new IntegratorMCFasterer(potentialMaster, random, params.temperature, box);

        MCMoveAtomFasterer mcMoveAtom = new MCMoveAtomFasterer(random, potentialMaster, box);
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

        LJMC3DFasterer sim = new LJMC3DFasterer(params);
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
        MeterPotentialEnergyFromIntegratorFasterer meterPE = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
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

            LJMC3DFasterer sim = new LJMC3DFasterer(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);
            MeterPotentialEnergyFromIntegratorFasterer meterPE = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
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
