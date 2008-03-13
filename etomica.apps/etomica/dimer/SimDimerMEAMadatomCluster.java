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
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";
    	
    	//Simulation 1 - MD and Dimer search
    	final SimDimerMEAMadatom sim1 = new SimDimerMEAMadatom(fileName, false, false, false, false, false, false);
    	sim1.activityIntegrateMD.setMaxSteps(mdSteps);
    	sim1.activityIntegrateDimer.setMaxSteps(2000);
    	
    	XYZWriter xyzwritermd = new XYZWriter(sim1.box);
        xyzwritermd.setFileName(fileName+"-md.xyz");
        xyzwritermd.setIsAppend(true);
    	XYZWriter xyzwriterd = new XYZWriter(sim1.box);
        xyzwriterd.setFileName(fileName+"-d.xyz");
        xyzwriterd.setIsAppend(true);
        sim1.integratorMD.addIntervalAction(xyzwritermd);
        sim1.integratorMD.setActionInterval(xyzwritermd, 10);
        sim1.integratorDimer.addIntervalAction(xyzwriterd);
        sim1.integratorDimer.setActionInterval(xyzwriterd, 1);
        
        sim1.getController().actionPerformed();
        /**
        //Simulation 2 - Fine grain Dimer search
        final SimDimerMEAMadatom sim2 = new SimDimerMEAMadatom(fileName, true, false, false, false);
        sim2.activityIntegrateDimer.setMaxSteps(100);
        sim2.getController().actionPerformed();
        
        //Simulation 3 - Vibrational normal mode analysis
        final SimDimerMEAMadatom sim3 = new SimDimerMEAMadatom(fileName+"_fine_saddle", false, true, false, false);
        sim3.getController().actionPerformed();
        
        //Simulation 4 - Minimum Search - A direction
        final SimDimerMEAMadatom sim4 = new SimDimerMEAMadatom(fileName, false, false, true, false);
    	sim4.activityIntegrateMin.setMaxSteps(500);
        sim4.getController().actionPerformed();
        
	    //Simulation 5 - Vibrational normal mode analysis
        final SimDimerMEAMadatom sim5 = new SimDimerMEAMadatom(fileName+"_A_minimum", false, true, false, false);
        sim5.getController().actionPerformed();
        
        //Simulation 6 - Minimum Search - B direction
        final SimDimerMEAMadatom sim6 = new SimDimerMEAMadatom(fileName, false, false, true, true);
    	sim6.activityIntegrateMin.setMaxSteps(500);
        sim6.getController().actionPerformed();
        
	    //Simulation 7 - Vibrational normal mode analysis
        final SimDimerMEAMadatom sim7 = new SimDimerMEAMadatom(fileName+"_B_minimum", false, true, false, false);
        sim7.getController().actionPerformed();
        */
    }

}
