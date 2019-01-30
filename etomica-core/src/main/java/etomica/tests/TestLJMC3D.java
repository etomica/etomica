/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
public class TestLJMC3D extends Simulation {

    public IntegratorMCFasterer integrator;
    public MCMoveAtomFasterer mcMoveAtom;
    public SpeciesSpheresMono species;
    public Box box;
    public PotentialMasterCellFasterer potentialMaster;
    public P2LennardJones potential;
    public Controller controller;

    public TestLJMC3D(int numAtoms, int numSteps, Configuration config) {
        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        double sigma = 1.0;
        box = this.makeBox();
        potentialMaster = new PotentialMasterCellFasterer(this, box, 2);
        integrator = new IntegratorMCFasterer(this, potentialMaster, box);
        mcMoveAtom = new MCMoveAtomFasterer(random, potentialMaster, box);
        mcMoveAtom.setStepSize(0.275 * sigma);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setTunable(false);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(numSteps);
        getController().addAction(activityIntegrate);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);

        double truncationRadius = 3.0 * sigma;

        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        AtomType leafType = species.getLeafType();
        potentialMaster.setPairPotential(leafType, leafType, potentialTruncated);

        config.initializeCoordinates(box);

        potentialMaster.init();
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3D.class);

        TestLJMC3D sim = new TestLJMC3D(numAtoms, params.numSteps, config);

        MeterPressureFasterer pMeter = new MeterPressureFasterer(sim.box, sim.potentialMaster);
        pMeter.setTemperature(sim.integrator.getTemperature());
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(10);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 2 * numAtoms);
        sim.integrator.getEventManager().addListener(pPump);
        MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 10);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().actionPerformed();
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);

        System.out.println("Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());

        double Z = ((DataDouble) ((DataGroup) pAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        avgPE /= numAtoms;
        System.out.println("Z=" + Z);
        System.out.println("PE/epsilon=" + avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k=" + Cv);

        if (Double.isNaN(Z) || Math.abs(Z + 0.25 - 0.4) > 0.2) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE + 4.56 - 0.2) > 0.04) {
            System.exit(2);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv - 0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(3);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 400000;
    }
}
