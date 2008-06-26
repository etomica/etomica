package etomica.models.oneDHardRods;

import etomica.normalmode.SimTarget;
import etomica.simulation.Simulation;
import etomica.space.Space;

public class Sim1DHR extends Simulation {


	private static final String APP_NAME = "Sim 1DHR";

	public Sim1DHR(Space _space){
        super(_space, true);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
        //set up simulation parameters
        int D = 1;
        int numMolecules = 10;
        double density = 0.5;
        double harmonicFudge = 1.0;
        double simTime = 1000;
        double temperature = 1;

        String filename = "normal_modes_1DHR_"+density+"_"+simTime+"_"+numMolecules+"_"+harmonicFudge;
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            numMolecules = Integer.parseInt(args[3]);
        }
        if (args.length > 4) {
            harmonicFudge = Double.parseDouble(args[4]);
        }
        
        System.out.println("Running 1D hard rod simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println(simTime+" time units");
        System.out.println("output data to "+filename);

        //instantiate simulation
        SimTarget sim = new SimTarget(Space.getInstance(D), numMolecules, density);
        
        
        
        
        
        
        
        System.out.println("Fini.");
	}

}
