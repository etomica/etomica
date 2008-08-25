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
	        
        String fileName = "bill";//args[0];
        
    	final String APP_NAME = "SimDimerMEAMGBCluster";
    	
    	final SimDimerMEAMGB sim = new SimDimerMEAMGB(new int[] {1,0,1});
    	
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
        sim.setMovableAtomsList();
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */
        
        //sim.initializeConfiguration(fileName+"_saddle");
        
        //sim.enableMolecularDynamics(1);
        
        sim.enableDimerSearch(fileName, 1500, false, false);
        
        //sim.enableMinimumSearch(fileName, false);
        
        /*
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"_saddle.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 5);
        */
        
        /*
        XYZWriter xyzwriterMin = new XYZWriter(sim.box);
        xyzwriterMin.setFileName(fileName+"_initial.xyz");
        xyzwriterMin.setIsAppend(true);
        sim.integratorMD.addIntervalAction(xyzwriterMin);
        sim.integratorMD.setActionInterval(xyzwriterMin, 1);
        */
        
        sim.getController().actionPerformed();

    }

}
