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
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
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
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    
    public TestLJMC3D(int numAtoms, int numSteps, Configuration config) {
        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);
        double sigma = 1.0;
        box = this.makeBox();
        integrator = new IntegratorMC(this, potentialMaster, box);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        mcMoveAtom.setStepSize(0.2 * sigma);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setTunable(false);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(numSteps);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);
        double truncationRadius = 3.0 * sigma;
        if (truncationRadius > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potentialTruncated.getRange());
        AtomType leafType = species.getLeafType();
        potentialMaster.addPotential(potentialTruncated, new AtomType[]{leafType, leafType});
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        config.initializeCoordinates(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();
//        WriteConfiguration writeConfig = new WriteConfiguration("LJMC3D"+Integer.toString(numAtoms),box,1);
//        integrator.addListener(writeConfig);
    }
 
    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3D.class);

        TestLJMC3D sim = new TestLJMC3D(numAtoms, params.numSteps, config);

        MeterPressure pMeter = new MeterPressure(sim.space);
        pMeter.setIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(10);
        DataPump pPump = new DataPump(pMeter,pAccumulator);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pPump);
        pumpListener.setInterval(2*numAtoms);
        sim.integrator.getEventManager().addListener(pumpListener);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyManager));
        
        sim.getController().actionPerformed();

        double Z = ((DataDouble) ((DataGroup) pAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        avgPE /= numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        
        if (Double.isNaN(Z) || Math.abs(Z+0.25) > 0.2) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+4.56) > 0.04) {
            System.exit(1);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv-0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(1);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 200000;
    }
}
