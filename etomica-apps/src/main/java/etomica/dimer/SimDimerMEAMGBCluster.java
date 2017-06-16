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

public class SimDimerMEAMGBCluster extends Simulation{
	

	public SimDimerMEAMGBCluster(Space space) {
		super(space);
	}

	public static void main(String[] args){
		
		String fileName = "test12"; //args[0];
        //int mdSteps = 10;//Integer.parseInt(args[1]);
        //int h = Integer.parseInt(args[1]);
        //int k = Integer.parseInt(args[2]);
        //int l = Integer.parseInt(args[3]);
        
        //int x = Integer.parseInt(args[4]);
        //int y = Integer.parseInt(args[5]);
        //int z = Integer.parseInt(args[6]);
        
        //double num = Double.parseDouble(args[1]);
        
        final String APP_NAME = "SimDimerMEAMGBCluster";
        
        final SimDimerMEAMGB sim = new SimDimerMEAMGB(new int[] {1,0,1}, new int[] {4,4,12});
        
        sim.initializeConfiguration("sngb101-4412-md");
        
        Vector dimerCenter = sim.getSpace().makeVector();
        dimerCenter.setX(0, sim.box.getBoundary().getBoxSize().getX(0)/2.0);
        dimerCenter.setX(1, 1.0);
        dimerCenter.setX(2, 0.0);
        Vector cubeSize = sim.getSpace().makeVector();
        cubeSize.setX(0, 6.0);
        cubeSize.setX(1, 8.0);
        cubeSize.setX(2, 8.0);
        
        if(sim.millerPlane[2] == 0){
            dimerCenter.setX(1, sim.box.getBoundary().getBoxSize().getX(1)/2.0);
            dimerCenter.setX(0, 1.0);
            dimerCenter.setX(2, 0.0);
            cubeSize.setX(0, 6.0);
            cubeSize.setX(1, 8.0);
            cubeSize.setX(2, 8.0);
        }

        
        
        sim.setMovableAtomsSphere(6.0, dimerCenter);
        //sim.setMovableAtomsCube(cubeSize, dimerCenter);
        //sim.setMovableAtomsList();
        
        //sim.removeAtoms(num, dimerCenter);
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */
        
        //sim.initializeConfiguration(fileName+"_saddle");
        
        //sim.enableMolecularDynamics(10000);
        
        sim.enableDimerSearch(fileName, 2500, false, false);
        sim.integratorDimer.setRotNum(1);
        
        
        //sim.enableMinimumSearch(fileName, true);
        
        /*
        //Limit MSD calculation to a specific species
        AtomIteratorFiltered aif = AtomIteratorFiltered.makeIterator(new AtomIteratorLeafAtoms(sim.box), 
        		new AtomFilterTypeInstance(sim.dimer.getChildType(0)));
        MeterMeanSquareDisplacement msd = new MeterMeanSquareDisplacement(sim.getSpace(), sim.integratorDimerMin);
        msd.setIterator((AtomIteratorBoxDependent)aif);
        
        
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"_B_minimum.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimerMin.addIntervalAction(xyzwriter);
        sim.integratorDimerMin.setActionInterval(xyzwriter, 5);
        */
        
        /*
        WriteConfiguration writer = new WriteConfiguration(sim.getSpace());
        writer.setBox(sim.box);
        writer.setConfName(fileName+"-md");
        sim.integratorMD.addIntervalAction(writer);
        sim.integratorMD.setActionInterval(writer, 10000);
		*/
           
        sim.getController().actionPerformed();
        
        /*
        Vector [] msdarray = msd.getDataAsArray();
        aif.reset();
        int i=0;
        System.out.println("-----MSD Info-----");
        for(IAtom a=aif.nextAtom(); a!=null; a=aif.nextAtom()){
        	System.out.println(((IAtomLeaf)a).getLeafIndex()+"     "+Math.sqrt(msdarray[i].squared()));
        	i++;
        }
		*/
    }

}
