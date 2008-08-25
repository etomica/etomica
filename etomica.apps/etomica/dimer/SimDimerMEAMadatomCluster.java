package etomica.dimer;

import etomica.action.WriteConfiguration;
import etomica.action.XYZWriter;
import etomica.api.IVector;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerMEAMadatomCluster extends Simulation{

	public SimDimerMEAMadatomCluster(Space space) {
		super(space);
	}

	public static void main(String[] args){
	        
        String fileName = args[0];
        //int mdSteps = Integer.parseInt(args[1]);
        //boolean ortho = Boolean.parseBoolean(args[2]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";

    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
    	IVector vect = sim.getSpace().makeVector();
        vect.setX(0, 10.0);
        vect.setX(1, 0.1);
        vect.setX(2, -1.0);
        
        sim.setMovableAtoms(80.0, vect);
        sim.setPotentialListAtoms();
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */
        
        sim.initializeConfiguration("sns-md-6rc-l");
        //sim.initializeConfiguration(fileName+"_saddle");
        
        sim.enableMolecularDynamics(1000);
        
        sim.enableDimerSearch(fileName, 1500, false, false);
        
        //sim.enableMinimumSearch(fileName, true);
        
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"_saddle.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 5);
        
        /*
        XYZWriter xyzwriterMin = new XYZWriter(sim.box);
        xyzwriterMin.setFileName(fileName+"_B_minimum.xyz");
        xyzwriterMin.setIsAppend(true);
        sim.integratorDimerMin.addIntervalAction(xyzwriterMin);
        sim.integratorDimerMin.setActionInterval(xyzwriterMin, 1);
        */
        
    	sim.getController().actionPerformed();

    }

}
