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
	        
        String fileName = args[0];
        
    	final String APP_NAME = "SimDimerMEAMGBCluster";
    	
    	final SimDimerMEAMGB sim = new SimDimerMEAMGB(fileName,new int[] {1,0,1});
    	
        IVector dimerCenter = sim.getSpace().makeVector();
        dimerCenter.setX(0, sim.box.getBoundary().getDimensions().x(0)/2.0);
        dimerCenter.setX(1, 1.0);
        dimerCenter.setX(2, 0.0);
        
        sim.initializeConfiguration("sn101-md");
        sim.setMovableAtoms(1.0, dimerCenter);
        dimerCenter.setX(1, 4.0);
        sim.removeAtoms(1.0, dimerCenter);
        sim.setMovableAtoms(2.5, dimerCenter);
        dimerCenter.setX(2, -2.0);
        sim.setMovableAtoms(2.0, dimerCenter);
        
        sim.initializeConfiguration("sn101-md11");
        
        sim.setMovableAtomsList();
        
        sim.enableMolecularDynamics(1000);
        
        sim.enableDimerSearch(fileName, 1000, true, false);
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+".xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 50);
        
        sim.getController().actionPerformed();

    }

}
