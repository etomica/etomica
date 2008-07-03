package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.SimTarget;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

public class Sim1DHR extends Simulation {

	int nA;
	double density;
	double temperature;
	String filename;
	double harmonicFudge;
	
	Integrator integrator;
	ActivityIntegrate activityIntegrate;
	IBox box;
	Boundary boundary;
	CoordinateDefinition coordinateDefinition;
	SpeciesSpheresMono species;
	
	
	NormalModes1DHR nm;
	private static final String APP_NAME = "Sim 1DHR";

	public Sim1DHR(Space _space, int numAtoms, double den, double tem, String file, double hF){
        super(_space, true);
        
        nA = numAtoms;
        density = den;
        temperature = tem;
        filename = file;
        harmonicFudge = hF; 
        
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        
        
        nm = new NormalModes1DHR(1);
        
        
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
