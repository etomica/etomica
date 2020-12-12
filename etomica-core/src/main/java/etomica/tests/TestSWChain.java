/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;


import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureHardFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHardFasterer;
import etomica.nbr.list.NeighborListManagerFasterer;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.util.ArrayList;
import java.util.List;

/**
 * Simple square-well chain simulation.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
 
public class TestSWChain extends Simulation {

    public IntegratorHardFasterer integrator;
    public Box box;
    static int chainLength = 10;

    public TestSWChain(Space _space, int numMolecules, Configuration config) {
        super(_space);
        setRandom(new RandomMersenneTwister(4));

        int numAtoms = numMolecules * chainLength;
        double sigma = 1.0;
        double sqwLambda = 1.5;
        double epsilon = 0.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        double timeStep = 0.005;
        double nbrRange = neighborRangeFac * sqwLambda * sigma;

        SpeciesGeneral species = new SpeciesBuilder(space)
                .withConformation(new ConformationLinear(space, sigma))
                .addCount(AtomType.simpleFromSim(this), chainLength)
                .setDynamic(true)
                .build();
        addSpecies(species);
        box = this.makeBox();

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(this);
        NeighborListManagerFasterer neighborManager = new NeighborListManagerFasterer(this, box, 2, nbrRange, bondingInfo);
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair potentialMaster = new PotentialComputePair(this, box, neighborManager);

        P2HardGeneric p2Bond = new P2HardGeneric(new double[]{sigma - bondFactor, sigma + bondFactor}, new double[]{Double.POSITIVE_INFINITY, 0, Double.POSITIVE_INFINITY});
        List<int[]> bondedIndices = new ArrayList<>();
        for (int i = 0; i < chainLength - 1; i++) {
            bondedIndices.add(new int[]{i, i + 1});
        }
        bondingInfo.setBondingPotentialPair(species, p2Bond, bondedIndices);

        // makes eta = 0.35
        double l = 14.4094 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        integrator = new IntegratorHardFasterer(potentialMaster, neighborManager, random, 0.01, 1.0, box, bondingInfo);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);

        P2HardGeneric potential = P2SquareWell.makeFastPotential(sigma, sqwLambda, epsilon);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), potential);

        box.getBoundary().setBoxSize(Vector.of(l, l, l));
        box.setNMolecules(species, numMolecules);
        config.initializeCoordinates(box);
        new BoxImposePbc(box, space).actionPerformed();
    }
    
    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numMolecules = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("SWChain%d.pos", numMolecules), TestSWChain.class);

        Space sp = Space3D.getInstance();
        TestSWChain sim = new TestSWChain(sp, numMolecules, config);
        int nSteps = (int) (params.numSteps / (numMolecules * chainLength * sim.integrator.getTimeStep()));
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps / 10));
        sim.integrator.resetStepCount();

        MeterPressureHardFasterer pMeter = new MeterPressureHardFasterer(sim.integrator);
        MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed();
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps));
        System.out.println("time: " + (System.nanoTime() - t1) / 1e9);

        double Z = pMeter.getDataAsScalar() * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x;
        avgPE /= numMolecules;
        System.out.println("Z=" + Z);
        System.out.println("PE/epsilon=" + avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numMolecules;
        System.out.println("Cv/k="+Cv);

        if (Double.isNaN(Z) || Math.abs(Z-4.5) > 1.5) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+19.32) > 0.12) {
            System.exit(1);
        }
        // actual value ~2
        if (Double.isNaN(Cv) || Cv < 0.5 || Cv > 4.5) {
            System.exit(1);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 100000;
    }
}
