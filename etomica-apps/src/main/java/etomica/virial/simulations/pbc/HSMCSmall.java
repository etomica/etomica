/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.pbc;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterWidomInsertion;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Hard sphere Monte Carlo Simulation of two species in 2D.
 * <p>
 * Species have the same sphere diameter but non-additive cross interactions.
 * (sigma12 = 1.2)
 *
 * @author David Kofke
 */

public class HSMCSmall extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public P2HardGeneric potential;
    public int chs;

    public HSMCSmall(int nAtoms, double L) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularPeriodic(getSpace(), L));
        PotentialMaster potentialMaster = new PotentialMaster(getSpeciesManager(), box, BondingInfo.noBonding());
        integrator = new IntegratorMC(potentialMaster, this.getRandom(), 1.0, box);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        box.setNMolecules(species, nAtoms);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        chs = 50;
        double sigma = chs * 0.01;
        potential = P2HardSphere.makePotential(sigma);
        AtomType type1 = species.getLeafType();
        potentialMaster.setPairPotential(type1, type1, potential);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

    }

    public Activity makeInitConfigActivity() {
        return new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
            boolean success = false;
            PotentialMaster potentialMaster = (PotentialMaster) integrator.getPotentialCompute();
            for (; chs <= 100; chs++) {
                potential.setCollisionDiameter(0, chs * 0.01);
                double u = potentialMaster.computeAll(false);
//                    System.out.println("chs "+chs*0.01+" "+u);
                if (u == Double.POSITIVE_INFINITY) {
                    chs--;
                    potential.setCollisionDiameter(0, chs * 0.01);
                    break;
                }
                if (chs == 100) {
                    success = true;
                    potentialMaster.init();
                }
            }
            if (!success) {
                long attempts = 0;
                while (chs < 100) {
                    attempts++;
                    integrator.reset();
                    for (int i = 0; i < 1000; i++) {
                        handle.yield(integrator::doStep);
                    }
                    chs++;
                    potential.setCollisionDiameter(0, chs * 0.01);
                    double u = potentialMaster.computeAll(false);
//                        System.out.println(attempts+"  chs "+chs*0.01+" "+u);
                    if (u == Double.POSITIVE_INFINITY) {
                        chs--;
                        potential.setCollisionDiameter(0, chs * 0.01);
                    }
                    if (attempts == 100000) {
                        throw new RuntimeException("Failed to find non-overlapping config");
                    }
                }
                potentialMaster.init();
            }
            integrator.reset();
            integrator.resetStepCount();
            }
        };
    }


    public static void main(String[] args) {
        SimParameters params = new SimParameters();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.L = 4;
            params.nAtoms = 5;
            params.steps = 10000000;
        }
        int nAtoms = params.nAtoms;
        double L = params.L;
        final HSMCSmall sim = new HSMCSmall(nAtoms, L);

        long steps = params.steps;

        System.out.println("HS MC");
        System.out.println("N: "+nAtoms);
        System.out.println("L: "+L);
        System.out.println("steps: "+steps);

        sim.getController().runActivityBlocking(sim.makeInitConfigActivity());

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps/10));

        int nBlocks = 100;
        int interval = 10;
        long blockSize = Math.max(steps / (nBlocks*interval), 1);

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(sim.box, sim.getRandom(), sim.integrator.getPotentialCompute(), 1.0);
        meterWidom.setSpecies(sim.species);
        meterWidom.setNInsert(1);
        AccumulatorAverageFixed accWidom = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom, 10);
        sim.integrator.getEventManager().addListener(pumpWidom);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        IData dataWidom = accWidom.getData();
        double avgWidom = dataWidom.getValue(accWidom.AVERAGE.index);
        double errWidom = dataWidom.getValue(accWidom.ERROR.index);
        double corWidom = dataWidom.getValue(accWidom.BLOCK_CORRELATION.index);

        System.out.println("\nexp(-mu/kT): "+avgWidom+"   err: "+errWidom+"   cor: "+corWidom);

        long t2 = System.nanoTime();

        System.out.println("time: "+(t2-t1)/1e9);
    }

    public static class SimParameters extends ParameterBase {
        public int nAtoms = 5;
        public long steps = 1000000;
        public double L = 2.2;
    }
}
