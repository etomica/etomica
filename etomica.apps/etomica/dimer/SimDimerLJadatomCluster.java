package etomica.dimer;

import etomica.simulation.Simulation;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerLJadatomCluster extends Simulation{
	
	 public static void main(String[] args){
	        
	        String fileName = "lj700"; //args[0];
	        int mdSteps = 800; //Integer.parseInt(args[1]);
	        
	    	final String APP_NAME = "SimDimerLJadatomCluster";
	    	/*
	    	//Simulation 1 - MD and Dimer search
	    	final SimDimerLJadatom sim1 = new SimDimerLJadatom(fileName, false, false, false, false);
	    	sim1.activityIntegrateMD.setMaxSteps(mdSteps);
	    	sim1.activityIntegrateDimer.setMaxSteps(700);
	        sim1.getController().actionPerformed();
	  
	        //Simulation 2 - Fine grain Dimer search
	        final SimDimerLJadatom sim2 = new SimDimerLJadatom(fileName, true, false, false, false);
	        sim2.activityIntegrateDimer.setMaxSteps(100);
	        sim2.getController().actionPerformed();
	        */
	        //Simulation 3 - Vibrational normal mode analysis
	        final SimDimerLJadatom sim3 = new SimDimerLJadatom(fileName+"_fine_saddle", false, true, false, false);
	        sim3.getController().actionPerformed();
	        
	        //Simulation 4 - Minimum Search - A direction
	        final SimDimerLJadatom sim4 = new SimDimerLJadatom(fileName, false, false, true, false);
	    	sim4.activityIntegrateMin.setMaxSteps(500);
	        sim4.getController().actionPerformed();
	        
		    //Simulation 5 - Vibrational normal mode analysis
	        final SimDimerLJadatom sim5 = new SimDimerLJadatom(fileName+"_A_minimum", false, true, false, false);
	        sim5.getController().actionPerformed();
	        
	        //Simulation 6 - Minimum Search - B direction
	        final SimDimerLJadatom sim6 = new SimDimerLJadatom(fileName, false, false, true, true);
	    	sim6.activityIntegrateMin.setMaxSteps(500);
	        sim6.getController().actionPerformed();
	        
		    //Simulation 7 - Vibrational normal mode analysis
	        final SimDimerLJadatom sim7 = new SimDimerLJadatom(fileName+"_B_minimum", false, true, false, false);
	        sim7.getController().actionPerformed();
	     
	    }




}
