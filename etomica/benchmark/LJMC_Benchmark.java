package etomica.benchmark;

import etomica.*;
import etomica.performance.*;
import etomica.units.*;
//import etomica.graphics.*;

public class LJMC_Benchmark {
    
      public static void main(String[] args) {
        
        // 0 = standard
        // 1 = performance
        // 2 = molecule layer
        final int VERSION = 2;
        final int nCycles = 1000;
        
        Default.BOX_SIZE=30.;
        Default.ATOM_SIZE=3.0;
        Default.TRUNCATE_POTENTIALS=false;
        
        PotentialLJ potential1 = null;
        Species species;
        
        switch(VERSION) {
            default:
            case 0:
            case 2:
               Simulation.instance = new Simulation(new Space3D());
               break;
            case 1:
                Simulation.instance = new Simulation(new SpaceP(3));
                break;
        }
        
        Simulation sim = Simulation.instance;
        sim.random.setSeed(0);        
        
        switch(VERSION) {
            default:
            case 0:
            case 1:
                species = new SpeciesSpheresMono();
                break;
            case 2:
                species = new SpeciesSpheres();
                break;
        }
                
        species.setNMolecules(64);
        
        switch(VERSION) {
            default:
            case 0:
               P2LennardJones potential = new P2LennardJones();
               potential.setSpecies(species, species);
               potential.setIterator(new ApiIntraspecies1A());
               break;
            case 1:  //performance
                potential1 = new PotentialLJ(Simulation.instance.hamiltonian.potential,species.getNMolecules());
                potential1.setSigma(Default.ATOM_SIZE);
                potential1.setEpsilon(Default.POTENTIAL_WELL);

                PotentialCalculationEnergySumPerformance energySum = new PotentialCalculationEnergySumPerformance();
                Simulation.instance.setEnergySum(energySum);
                break;
            case 2:  //molecule layer
                Potential2Group group = new Potential2Group();
                P2LennardJones potential2 = new P2LennardJones(group);
                group.setIterator(new ApiIntraspecies1A());
                group.setSpecies(species,species);
                potential2.setIterator(new AtomPairIteratorSynthetic(Simulation.instance.space,
                    new AtomIteratorChildren(), new AtomIteratorChildren()));
                break;
        }       
        
        Phase phase = new Phase(sim);
        
        Controller controller = new Controller();
        IntegratorMC integrator = new IntegratorMC();
        integrator.setDoSleep(false);
	    integrator.setInterval(species.getNMolecules());
	    
//	    MCMoveMolecule mcmovemolecule = new MCMoveMolecule(integrator);
	    MCMoveAtom mcmovemolecule = new MCMoveAtom(integrator);
       
///	    MeterCycles meter2  = new MeterCycles();
	    MeterPotentialEnergy meter4 = new MeterPotentialEnergy();
	    meter4.setUpdateInterval(10);
	    
///	    meter2.setName("cycles");
///	    meter4.setName("potentialenergy");
	    
///	    MeterValuesToLogFile displayLog = new MeterValuesToLogFile("vrB.BE9.0.C7.txt");
///	    displayLog.setUpdateInterval(10/*10000 DAK*/);
///	    displayLog.setResetVal(20/*200000 DAK*/);
///	    displayLog.setReset(true);
	    
        //MeterValuesBenchMark displayLog = new MeterValuesBenchMark("VanRoijBiaBetaE2.0NonGUI.txt");
///	    integrator.addIntervalListener(displayLog);
	    
	    sim.elementCoordinator.go();
        
        if(VERSION == 1) potential1.setArray(phase);
        
        integrator.setMaxSteps(nCycles*species.getNMolecules());
   //     ((Controller)sim.controller(0)).start();

        System.out.println("AtomCount: "+phase.atomCount());
        System.out.println("Molecules: "+phase.moleculeCount());
        System.out.println("Starting");
        integrator.initialize();
        Stopwatch timer = new Stopwatch().start();
        integrator.run();
   //     integrator.start();
        timer.stop();
        System.out.println("Done");
        System.out.println(meter4.average());
        System.out.println("Elapsed time (sec): "+0.001*timer.getElapsedTime());
    }
}

/*
 
 RESULTS 
 
 2D average energy: -28125.923188136127
 3D average energy: - 9783.048798461032
 
 Time -- Configuration
 -----------------------
 26.9 -- DAK Thinkpad; 64 atoms; 1000 cycles; 3D (floor/ceil for nearestImage)
 16.0 -- DAK/64/1000; 2D (while-loops for nearest image); 
 18.8 -- same, 3D
 14.5 -- DAK/64/1000; 2D/while; pass action to iterator (must comment setFirst in AtomPairIterator)
 16.4 -- same, 3D
  8.3 -- DAK/64/1000; 3D/while; performance Space, PotentialLJ (avg energy = -9341.847638448702)
 15.8 -- DAK/64/1000; 2D/while; revised Boundary to store 0.5*dimensions
 16.8 -- same, 3D
 13.8 -- DAK/64/1000; 2D/while; removed momentum calculations in CoordinatePair.reset
 15.1 -- same, 3D
 14.8 -- DAK/64/1000; (0) 3D/while/no momentum; r2 computed as needed in CoordinatePair
 12.6 -- DAK/64/1000; 3D/while/nomom/r2asNeeded; no truncation in Potential2SoftSpherical
  5.3 -- DAK/64/1000; (1) 3D/performance; PBC in PotentialLJ (energy = -9341.8...702)
  6.4 -- DAK/64/1000; 3D/performance/ PBC in boundary, local boundary field, dr passed as 3D.Vector
  
  118 -- DAK/64/1000; (2) 3D/momentum; molecule layer (avg energy = -9583.830204672915)
 90.5 -- DAK/64/1000; 3D/noMom/molecule
 78.9 -- DAK/64/1000; 3D/mom/molecule; store depth info
 93.0 -- DAK/64/1000; 3D/mom/molecule; AtomTreeNodeRecursive
116.0 -- AtomTreeNodeRecursive
105
 88   -- first/last atom not in space
 74   -- store childrenAreGroups
 83   -- index in node
 38   -- (2) AtomIteratorChildren in pair iterator for potential2
 
 43.9 -- Home W2K Athlon/64/1000 (2)  (32.6 when run with HotSpot Server)
  2.7 -- Home/64/1000 (1)            (3.0 HotSpot Server, based on 10K cycles) (34.0 HotSpot Classic)
  6.7 -- Home/64/1000 (0)            (6.2 HotSpot Server, based on 20K cycles)
 */