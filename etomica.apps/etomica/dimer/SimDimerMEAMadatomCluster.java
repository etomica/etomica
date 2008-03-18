package etomica.dimer;

import etomica.action.XYZWriter;
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
        int mdSteps = Integer.parseInt(args[1]);
        boolean ortho = Boolean.parseBoolean(args[2]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";

    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
    	sim.initializeConfiguration(fileName);
    	sim.enableDimerSearch(fileName, mdSteps, ortho, false);

    	XYZWriter xyzwriter = new XYZWriter(sim.box);
    	xyzwriter.setFileName(fileName+".xyz");
    	xyzwriter.setIsAppend(true);
    	
    	sim.integratorDimer.addIntervalAction(xyzwriter);
    	sim.integratorDimer.setActionInterval(xyzwriter, 5);
    	
    	sim.getController().actionPerformed();

    }

}
