package etomica.simulations;

import etomica.*;

import etomica.units.*;

import etomica.datamanager.*;

import etomica.performance.*;

import etomica.graphics.*;

public class Lennard3DPerfGraphic{

      public static void main(String[] args) {
      
        SimulationGraphic sim = new SimulationGraphic(new SpaceP(3));

        // Default.IS_GRAPHIC=true;

        Default.BOX_SIZE=50.0;

        Default.ATOM_SIZE=3.0;

        Default.TRUNCATE_POTENTIALS=false;

        SimulationGraphic.instance = sim;

        Phase phase = new Phase(sim);

        

        SpeciesSpheresMono species = new SpeciesSpheresMono();

        species.setNMolecules(64);

        //SHOULD I DO THIS ?

	//YES

        PotentialLJ potential = new PotentialLJ(Simulation.instance.hamiltonian.potential,species.getNMolecules());

        
        
        potential.setSigma(Default.ATOM_SIZE);
        potential.setEpsilon(Default.POTENTIAL_WELL);

        PotentialCalculationEnergySumPerformance energySum = new PotentialCalculationEnergySumPerformance();
        
        sim.setEnergySum(energySum);
        
        

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

        //DisplayPhase displayphase = new DisplayPhase();
	    DisplayBox box1 = new DisplayBox();
	    DisplayBox box2 = new DisplayBox();
	    //DisplayBox box3 = new DisplayBox();
	    DisplayBox box4 = new DisplayBox();
	    
	    box1.setMeter(meter1);
	    box2.setMeter(meter2);
	    box2.setPrecision(10);
        box4.setMeter(meter4);
	    box4.setPrecision(10);
	    //box3.setUseCurrentValue(false);
	    box4.setUseCurrentValue(false);
	    
	    MeterValuesToLogFile displayLog = new MeterValuesToLogFile("ljpT.txt");
	    displayLog.setReset(true);
	    displayLog.setUpdateInterval(1000);
	    displayLog.setResetVal(200000);
	    displayLog.setCycleCount(1000);
	   
	    displayLog.setStartTime(System.currentTimeMillis());
        integrator.addIntervalListener(displayLog);
	    sim.elementCoordinator.go();

        potential.setArray(phase);
        SimulationGraphic.makeAndDisplayFrame(sim);
        ((Controller)sim.controller(0)).start();

        

    }

    

    

}

