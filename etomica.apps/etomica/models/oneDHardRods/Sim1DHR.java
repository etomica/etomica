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
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

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

        
        
        
        
        
        
        
        
        
        
        
        
       

        
        
        nm = new NormalModes1DHR(1);
        
        
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
		/*
		 * This whole setup defines a set of default parameters
		 * in the inner class Sim1DHRParams.  These parameters can be changed
		 * individually in an appropriately named file, without affecting
		 * the values of the other parameters.  The order of definition in the
		 * file is irrelevant.
		 * 
		 */
		Sim1DHRParams params = new Sim1DHRParams();
		String inputFilename = null;
		if(args.length > 0){
			inputFilename = args[0];
		}
        if(inputFilename != null){
        	ReadParameters readParameters = new ReadParameters(inputFilename, params);
        	readParameters.readParameters();
        }
        double density = params.density;
        double numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if(filename.length() ==0){
        	filename = "normal_modes_1DHR _" + numMolecules;
        }
        String refFileName = args.length>0 ? filename+"_ref" : null;
        
        System.out.println("Running 1D hard rod simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+ " total steps of 1000");
        System.out.println("output data to "+filename);
        
        
        
        
        //instantiate simulation
        SimTarget sim = new SimTarget(Space.getInstance(D), numMolecules, density);
        
        
        
        
        
        
        
        System.out.println("Fini.");
	}

	
    /**
     * Inner class for parameters understood by the Sim1DHR constructor
     */
    public static class Sim1DHRParams extends ParameterBase {
        public int numMolecules = 32;
        public double density = 0.5;
        public int D = 1;
        public long numSteps = 100000;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
    }
}
