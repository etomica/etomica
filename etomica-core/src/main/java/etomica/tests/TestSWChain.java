/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
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
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheres;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple square-well chain simulation.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
 
public class TestSWChain extends Simulation {

    public IntegratorHard integrator;
    public Box box;

    public TestSWChain(Space _space, int numMolecules, double simTime, Configuration config) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        int chainLength = 10;
        int numAtoms = numMolecules * chainLength;
        double sigma = 1.0;
        double sqwLambda = 1.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        double timeStep = 0.005;
        simTime /= chainLength;
        int nSteps = (int) (simTime / timeStep);

        // makes eta = 0.35
        double l = 14.4094 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        box = new Box(space);
        integrator = new IntegratorHard(this, potentialMaster, space, box);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(nSteps);
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(neighborRangeFac * sqwLambda * sigma);

        SpeciesSpheres species = new SpeciesSpheres(this, _space, chainLength);
        species.setIsDynamic(true);
        addSpecies(species);
        P2HardBond bonded = new P2HardBond(space, sigma, bondFactor, false);
        PotentialGroup potentialChainIntra = potentialMaster.makePotentialGroup(1);
        potentialChainIntra.addPotential(bonded, ApiBuilder.makeAdjacentPairIterator());

        potentialMaster.addPotential(potentialChainIntra, new ISpecies[]{species});
        ((ConformationLinear) species.getConformation()).setBondLength(sigma);

        P2SquareWell potential = new P2SquareWell(space, sigma, sqwLambda, 0.5, false);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});
        CriterionInterMolecular sqwCriterion = (CriterionInterMolecular) potentialMaster.getCriterion(potential);
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);
        addBox(box);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{l, l, l}));
        box.setNMolecules(species, numMolecules);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        integrator.setBox(box);
        config.initializeCoordinates(box);
    }
    
    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numMolecules = params.numAtoms;
        double simTime = params.numSteps/numMolecules;
        ConfigurationFile config = new ConfigurationFile("SWChain"+Integer.toString(numMolecules));

        Space sp = Space3D.getInstance();
        TestSWChain sim = new TestSWChain(sp, numMolecules, simTime, config);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed();
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyPump);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.box.getBoundary().volume()/(sim.box.getMoleculeList().getMoleculeCount()*sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        avgPE /= numMolecules;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numMolecules;
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
