/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class LJMCNVT extends Simulation {

    public final SpeciesSpheresMono species;
    public final Box box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtom mcMoveAtom;

    public LJMCNVT(int numMolceules, double density, double temperature, long numSteps) {
        super(Space3D.getInstance());
        int N = numMolceules;  //number of atoms

        System.out.println("numAtom=" +numMolceules);
        System.out.println("temperature = "+temperature);
        System.out.println("numSteps = "+numSteps);
        System.out.println("initial density = "+density+"\n");

        //controller and integrator
	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
	    species.setIsDynamic(true);
        addSpecies(species);

        //construct box
	    box = new Box(space);
        addBox(box);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);


        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, 4, space);
        potentialMaster.setCellRange(2);
        integrator = new IntegratorMC(potentialMaster, random, temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(numSteps);
        getController().addAction(activityIntegrate);

        //instantiate several potentials for selection in combo-box
        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential,2.5);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();

        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMinAdjustStep(1.05);

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) {

        LJMCParam params = new LJMCParam();
        ParseArgs.doParseArgs(params, args);
        int numMolecules = params.numMolecules;
        double density = params.density;
        double temperature = params.temperature;
        long numSteps = params.numSteps;
            
        LJMCNVT sim = new LJMCNVT(numMolecules,density,temperature,numSteps);
        sim.getController().actionPerformed();

        long t1 = System.currentTimeMillis();

        final MeterPotentialEnergyFromIntegrator meterE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageCovariance energyAccumulator = new AccumulatorAverageCovariance(10);
        DataPumpListener energyManager = new DataPumpListener(meterE, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyManager);

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().reset();
        sim.getController().actionPerformed();

        double avgEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index).getValue(0);
        System.out.println("average enthalpy (sim) = "+avgEnthalpy);

        double varEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(0);

        System.out.println("varEnthalpy = "+varEnthalpy );

        double avgCv = varEnthalpy/(temperature*temperature*numMolecules);

        System.out.println("average Cv(sim) = "+avgCv+"\n");

        long t2 = System.currentTimeMillis();
        System.out.println("time = "+(t2-t1)/1000.0+"\n");
    }

    public static class LJMCParam extends ParameterBase {
        public int numMolecules = 256;
        public double density = 0.082226;//1g/cm3=1000/18.02mol/L
        public double temperature = 1.5;//Kelvin
        public long numSteps = 1000000;
    }

}
