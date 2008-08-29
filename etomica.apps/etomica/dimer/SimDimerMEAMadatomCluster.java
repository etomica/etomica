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
        boolean ortho = Boolean.parseBoolean(args[1]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";

    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
    	IVector vect = sim.getSpace().makeVector();
        vect.setX(0, 10.0);
        vect.setX(1, 0.1);
        vect.setX(2, -1.0);
        
        sim.setMovableAtoms(90.0, vect);
        sim.setPotentialListAtoms();
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */
        
        sim.initializeConfiguration("sns-initial");
        //sim.initializeConfiguration(fileName+"_saddle");
        
        sim.enableMolecularDynamics(1000);
        
        sim.enableDimerSearch(fileName, 1500, ortho, false);
        
        //sim.enableMinimumSearch(fileName, true);
        /*
        WriteConfiguration poswriter = new WriteConfiguration(sim.getSpace());
        poswriter.setConfName(fileName);
        poswriter.setBox(sim.box);
        sim.integratorDimer.addIntervalAction(poswriter);
        sim.integratorDimer.setActionInterval(poswriter, 1000);
        */
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+".xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 1);
        
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
