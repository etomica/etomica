package etomica.simulations;

import etomica.*;

import etomica.units.*;

//import etomica.datamanager.*;

import etomica.performance.*;

//import etomica.graphics.*;

public class Lennard3DPerformance{

      public static void main(String[] args) {

        Simulation sim = new Simulation(new SpaceP(3));

       

        Default.BOX_SIZE=50;

        Default.ATOM_SIZE=3.0;

        Default.TRUNCATE_POTENTIALS=false;

        Simulation.instance = sim;

        Phase phase = new Phase(sim);

        

        SpeciesSpheresMono species = new SpeciesSpheresMono();

        species.setNMolecules(64);

        //SHOULD I DO THIS ?

	//YES

        PotentialBase potential = new PotentialBase(Simulation.instance.hamiltonian.potential,species.getNMolecules());

        
        
        potential.setSigma(Default.ATOM_SIZE);
        potential.setEpsilon(Default.POTENTIAL_WELL);

///        PotentialCalculationEnergySumPerformance energySum = new PotentialCalculationEnergySumPerformance();
        
///        sim.setEnergySum(energySum);
        
        

        Controller controller = new Controller();

        IntegratorMC integrator = new IntegratorMC();

        integrator.setDoSleep(false);

	    integrator.setInterval(species.getNMolecules());

	    

	    MCMoveAtom mcmoveAtom = new MCMoveAtom(integrator);

        

        MeterDensity meter1 = new MeterDensity();

	    MeterCycles meter2  = new MeterCycles();

	    MeterPotentialEnergy meter4 = new MeterPotentialEnergy();

	    

	    meter1.setName("density1");

	    meter2.setName("cycles");

	    meter4.setName("potentialenergy");

///	    MeterValuesToLogFile displayLog = new MeterValuesToLogFile("ljpT.txt");
/*	    displayLog.setReset(true);
	    displayLog.setUpdateInterval(100);
	    displayLog.setResetVal(200000);
	    displayLog.setCycleCount(100);
	   
	    displayLog.setStartTime(System.currentTimeMillis());
      integrator.addIntervalListener(displayLog);
  */	    sim.elementCoordinator.go();

        potential.setArray(phase);

        ((Controller)sim.controller(0)).start();

        

    }

    

    

}

