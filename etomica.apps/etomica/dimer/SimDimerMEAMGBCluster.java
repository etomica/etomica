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

public class SimDimerMEAMGBCluster extends Simulation{
	

	public SimDimerMEAMGBCluster(Space space) {
		super(space);
	}

	public static void main(String[] args){
	        
        String fileName = "sn201";//args[0];
        //int mdSteps = Integer.parseInt(args[1]);
        //int h = Integer.parseInt(args[1]);
        //int k = Integer.parseInt(args[2]);
        //int l = Integer.parseInt(args[3]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";
    	
    	final SimDimerMEAMGB sim = new SimDimerMEAMGB(fileName,new int[] {2,0,1});
    	
        IVector dimerCenter = sim.getSpace().makeVector();
        dimerCenter.setX(0, sim.box.getBoundary().getDimensions().x(0)/2.0);
        dimerCenter.setX(1, 1.0);
        dimerCenter.setX(2, 0.0);
    	    	
    	sim.initializeConfiguration(fileName);
    	sim.setMovableAtoms(5.0, dimerCenter);
    	sim.setMovableAtomsList();
    	sim.removeAtoms(2.0, dimerCenter);
    	sim.enableDimerSearch(fileName, 1000, false, false);
    	
    	XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"-dimer.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 1);
        
        WriteConfiguration config = new WriteConfiguration(sim.getSpace());
        config.setBox(sim.box);
        config.setConfName(fileName+"-dimer");
        sim.integratorDimer.addIntervalAction(config);
        sim.integratorDimer.setActionInterval(config,10);
        
        sim.getController().actionPerformed();

    }

}
