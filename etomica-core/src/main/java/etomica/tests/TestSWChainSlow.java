/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBond;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple square-well chain simulation.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */

public class TestSWChainSlow extends Simulation {

    public IntegratorHard integrator;
    public Box box;
    static int chainLength = 10;

    public TestSWChainSlow(Space _space, int numMolecules, Configuration config) {
        super(_space);

        SpeciesGeneral species = new SpeciesBuilder(space)
                .withConformation(new ConformationLinear(space))
                .addCount(AtomType.simpleFromSim(this), chainLength)
                .setDynamic(true)
                .build();
        addSpecies(species);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        int numAtoms = numMolecules * chainLength;
        double sigma = 1.0;
        double sqwLambda = 1.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        double timeStep = 0.005;

        // makes eta = 0.35
        double l = 14.4094 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(neighborRangeFac * sqwLambda * sigma);
        P2HardBond bonded = new P2HardBond(space, sigma, bondFactor, false);
        PotentialGroup potentialChainIntra = potentialMaster.makePotentialGroup(1);
        potentialChainIntra.addPotential(bonded, ApiBuilder.makeAdjacentPairIterator());

        potentialMaster.addPotential(potentialChainIntra, new ISpecies[]{species});
        ((ConformationLinear) species.getConformation()).setBondLength(sigma);

        P2SquareWell potential = new P2SquareWell(space, sigma, sqwLambda, 0.5, false);

        AtomType sphereType = species.getAtomType(0);
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});
        CriterionInterMolecular sqwCriterion = (CriterionInterMolecular) potentialMaster.getCriterion(sphereType, sphereType);
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);
        box.getBoundary().setBoxSize(Vector.of(l, l, l));
        box.setNMolecules(species, numMolecules);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numMolecules = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("SWChain%d.pos", numMolecules), TestSWChainSlow.class);

        Space sp = Space3D.getInstance();
        TestSWChainSlow sim = new TestSWChainSlow(sp, numMolecules, config);
        int nSteps = (int) (params.numSteps / (numMolecules * chainLength * sim.integrator.getTimeStep()));
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps / 10));
        sim.integrator.resetStepCount();

        long bs = nSteps / 100;
        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage uAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, uAccumulator);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps));
        System.out.println("time: " + (System.nanoTime() - t1) / 1e9);

        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + avgP + " " + errP + " " + corP);

        double avgU = uAccumulator.getData(pAccumulator.AVERAGE).getValue(0) / numMolecules;
        double errU = uAccumulator.getData(pAccumulator.ERROR).getValue(0) / numMolecules;
        double corU = uAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("U " + avgU + " " + errU + " " + corU);

        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) uAccumulator.getData()).getData(uAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numMolecules;
        System.out.println("Cv/k " + Cv);

        // expected values based on 2x10^7 steps
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 3.008496525057149e-01 - 1.570490504414807e-01 / numMolecules;
        double stdevP = 0.006;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedU = -1.930034804482143e+01 + 1.424254285711056e+00 / numMolecules;
        double stdevU = 0.06;   // short sim stdev is smaller, but short avg deviates from long avg
        if (Double.isNaN(avgU) || Math.abs(avgU - expectedU) / stdevU > 4) {
            System.exit(2);
        }

        // Cv is 2.0, but can be pretty far off, especially for N=4000
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 200000;
    }
}
