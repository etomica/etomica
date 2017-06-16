/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.space.Vector;
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
        //boolean ortho = Boolean.parseBoolean(args[1]);
        
    	final String APP_NAME = "SimDimerMEAMadatomCluster";

    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
        Vector vect = sim.getSpace().makeVector();
        vect.setX(0, 9.8);
        vect.setX(1, -0.2);
        vect.setX(2, -0.2);
        
        sim.setMovableAtoms(100.0, vect);
        
        sim.setPotentialListAtoms();
        
        //sim.removeAtoms(2.9, vect);
        
        //sim.initializeConfiguration("sns101-initial3");
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */
        
        
        //sim.initializeConfiguration(fileName+"_saddle");
        
        //sim.enableMolecularDynamics(5000);
        
        sim.enableDimerSearch(fileName, 2000, false, false);
        sim.integratorDimer.setRotNum(0);
        
        //sim.enableMinimumSearch(fileName, true);
        
        /*
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"_B_minimum.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimerMin.addIntervalAction(xyzwriter);
        sim.integratorDimerMin.setActionInterval(xyzwriter, 2);
        */
        
    	sim.getController().actionPerformed();

    }

}
