package etomica.dimer;

import etomica.action.WriteConfiguration;
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

public class SimDimerMEAMGBCluster extends Simulation{
	

	public SimDimerMEAMGBCluster(Space space) {
		super(space);
	}

	public static void main(String[] args){
	        
        String fileName = args[0];
        int mdSteps = Integer.parseInt(args[1]);
        int h = Integer.parseInt(args[2]);
        int k = Integer.parseInt(args[3]);
        int l = Integer.parseInt(args[4]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";
    	
    	final SimDimerMEAMadatomGB sim = new SimDimerMEAMadatomGB("fileName",new int[] {h,k,l});
    	
    	XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+".xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorMD.addIntervalAction(xyzwriter);
        sim.integratorMD.setActionInterval(xyzwriter, 10);
        
        WriteConfiguration config = new WriteConfiguration(sim.getSpace());
        config.setBox(sim.box);
        config.setConfName(fileName);
        sim.integratorMD.addIntervalAction(config);
        sim.integratorMD.setActionInterval(config,mdSteps);
    
    	sim.activityIntegrateMD.setMaxSteps(mdSteps);	
    	sim.activityIntegrateDimer.setMaxSteps(0);
        sim.getController().actionPerformed();

    }

}
