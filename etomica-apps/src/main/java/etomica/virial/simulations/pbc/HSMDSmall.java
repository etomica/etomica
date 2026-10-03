/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.pbc;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterWidomInsertion;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2HardGeneric;
import etomica.potential.compute.NeighborManagerHard;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class HSMDSmall extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    /**
     * The Box holding the atoms.
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesGeneral species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardGeneric potential;

    public final PotentialComputePair potentialMaster;

    public HSMDSmall(int nAtoms, double L) {

        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularPeriodic(getSpace(), L));

        double sigma = 1.0;
        NeighborManagerHard neighborManager = new NeighborManagerSimpleHard(box);
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        potential = new P2HardGeneric(new double[]{sigma}, new double[]{Double.POSITIVE_INFINITY}, false);
        AtomType leafType = species.getLeafType();

        potentialMaster.setPairPotential(leafType, leafType, potential);

        integrator = new IntegratorHard(potentialMaster.getPairPotentials(), neighborManager, random, 0.01, 1, box, getSpeciesManager());
        integrator.setIsothermal(false);
        integrator.setMaxCollisionDiameter(species.getLeafType(), 1);

        box.setNMolecules(species, nAtoms);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
    }

    @Override
    public IntegratorHard getIntegrator() {
        return integrator;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD3DParam params = new HSMD3DParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.L = 2;
            params.nAtoms = 2;
            params.steps = 10000000;
        }
        final HSMDSmall sim = new HSMDSmall(params.nAtoms, params.L);

        if (true) {
            final String APP_NAME = "HSMD3D";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

            simGraphic.makeAndDisplayFrame(APP_NAME);
            return;
        }

        System.out.println("HS MD");
        System.out.println("N: "+params.nAtoms);
        System.out.println("L: "+params.L);
        System.out.println("steps: "+params.steps);

        long steps = params.steps;
        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps/10));
        sim.integrator.resetStepCount();

        int nBlocks = 100;
        int interval = 20;
        long blockSize = Math.max(steps / (nBlocks*interval), 1);

        MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, interval);
        sim.integrator.getEventManager().addListener(pumpP);

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(sim.box, sim.getRandom(), sim.potentialMaster, 1.0);
        meterWidom.setSpecies(sim.species);
        meterWidom.setNInsert(10);
        AccumulatorAverageFixed accWidom = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom, interval);
        sim.integrator.getEventManager().addListener(pumpWidom);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        IData dataP = accP.getData();
        double avgP = dataP.getValue(accP.AVERAGE.index);
        double errP = dataP.getValue(accP.ERROR.index);
        double corP = dataP.getValue(accP.BLOCK_CORRELATION.index);

        IData dataWidom = accWidom.getData();
        double avgWidom = dataWidom.getValue(accWidom.AVERAGE.index);
        double errWidom = dataWidom.getValue(accWidom.ERROR.index);
        double corWidom = dataWidom.getValue(accWidom.BLOCK_CORRELATION.index);

        System.out.println("P: "+avgP+"   err: "+errP+"   cor: "+corP);
        System.out.println("exp(-mu/kT): "+avgWidom+"   err: "+errWidom+"   cor: "+corWidom);
        long t2 = System.nanoTime();

        System.out.println("time: "+(t2-t1)/1e9);
    }

    public static HSMD3DParam getParameters() {
        return new HSMD3DParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class HSMD3DParam extends ParameterBase {
        public int nAtoms = 5;
        public double L = 4;
        public long steps = 1000000;
    }
}
