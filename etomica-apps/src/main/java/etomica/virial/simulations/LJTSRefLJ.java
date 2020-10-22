/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;





import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationClusterMove;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * The direct-sampling formulation of Mayer-Sampling Monte Carlo is employed to compute 
 * a virial coefficient of truncated-and-shifted Lennard Jones (LJTS), at a particular temperature.
 * 
 * Lennard-Jones is employed as the reference system.  
 * 
 * Virial coefficients of multiple LJTS can be evaluated simultaneously.  
 * 
 * Kate Shaul, 2009
 */
public class LJTSRefLJ {


    public static void main(String[] args) {

    	VirialLJParam params = new VirialLJParam();

        double temperature; final int nPoints; int log10Steps;
        
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperature = params.temperature;
            log10Steps = params.log10Steps;
            
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 3) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperature = Double.parseDouble(args[1]);
            log10Steps = Integer.parseInt(args[2]);
           
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments.");
        }

 
        long steps = (long)Math.pow(10, log10Steps);

        
        double ref = 0;
        double refError = 0;
        try {
        	//FileReader fileReader = new FileReader("/usr/users/andrew/virial/LJ/comp/B" + nPoints + ".dat");
        	FileReader fileReader = new FileReader("/san/user/schadel2/u2/andrewLJData/B" + nPoints + ".dat");
            BufferedReader bufReader = new BufferedReader(fileReader);
    
            String line = bufReader.readLine();
            while (line != null) {
            
                String[] valueStrings = line.split(" +");
                
                if (Double.parseDouble(valueStrings[0]) == temperature) {
                	ref = Double.parseDouble(valueStrings[1]);
                	
                	if (nPoints > 3) {
                		refError = Double.parseDouble(valueStrings[2]);
                	}
                	
                }
               
                line = bufReader.readLine();
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file ");
        }
        
        // Multiple LJTS (with different truncation distances, rc) can be considered by the same simulation
        double[] rc = new double[] {1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
        //double[] rc = new double[] {2.5};
        
        System.out.println("LJTS B" + nPoints);
        System.out.println("T* = " + temperature );
        if (rc.length == 1) {
        	System.out.println("rc = " + rc[0] );
        }
        System.out.println();
        System.out.println("Direct sampling with LJ as reference");
        System.out.println("LJ B" + nPoints + " at T* = " + ref + " +/- " + refError);
        System.out.println();
        System.out.println("10^" + log10Steps + " steps");
        System.out.println();
		
        Space space = Space3D.getInstance();
        
        
        Potential2Spherical pLJ = new P2LennardJones(space,1.0,1.0);
        
        MayerGeneralSpherical fRef = new MayerGeneralSpherical(pLJ);
        MayerESpherical eRef = new MayerESpherical(pLJ);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        
        
        ClusterAbstract[] targetCluster = new ClusterAbstract[rc.length];
        for (int i=0;i<rc.length;i++) {
        	
        	P2SoftSphericalTruncatedShifted pLJTS = new P2SoftSphericalTruncatedShifted(space, (Potential2SoftSpherical)pLJ, rc[i]);
        	
        	MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pLJTS);
            MayerESpherical eTarget = new MayerESpherical(pLJTS);
            
            targetCluster[i]= Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
            targetCluster[i].setTemperature(temperature);
            
        }
        
        
        
        //////////////////////////////////////////////////////////////////////////////////////////////
        // The sampleCluster is the system used to sample configuration space.  This should be the full LJ potential, as its important configuration space 
        // encompasses that of LJTS.
        //////////////////////////////////////////////////////////////////////////////////////////////
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        
        SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
		
        final SimulationVirial sim = new SimulationVirial(space, speciesFactory, temperature, sampleCluster,refCluster,targetCluster);
        
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box);
        
        
        
        sim.equilibrate(steps/40);
ActivityIntegrate ai = new ActivityIntegrate(sim.integrator, steps);
System.out.println("equilibration finished");

        System.out.println();
        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
sim.getController().runActivityBlocking(ai);

//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        
       
        
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        
        System.out.println();
        System.out.println("LJ average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
                           +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
                           +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);
        
        double[] ratio = new double[rc.length];
        double[] error = new double[rc.length];
        
        System.out.println();
        System.out.println("LJTS");
        System.out.println();
        
    	for (int i=0;i<rc.length;i++) {
    		
    		ratio[i] = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[i+1];
            error[i] = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[i+1];
        
    		System.out.println("rc = "+ rc[i] + " average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[i+1]
                               +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[i+1]
                               +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[i+1]);
       
        }
    	
    	System.out.println();
        
        for (int i=0;i<rc.length;i++) {
        	        
        	System.out.println("rc = "+ rc[i] + " ratio average: "+ratio[i]+", error: "+error[i]);
        	
        }
        
        System.out.println();
        
        if (nPoints > 3) {
        
        	for (int i=0;i<rc.length;i++) {
        	
        		System.out.println("rc = "+ rc[i] + " Bn: "+ratio[i]*ref+", error: "+ Math.sqrt(error[i]*error[i]*ref*ref + refError*refError*ratio[i]*ratio[i]));
        		
        	}
        		
    	} else {
    		
    		for (int i=0;i<rc.length;i++) {
        		
        		System.out.println("rc = "+ rc[i] + " Bn: "+ratio[i]*ref+", error: "+ error[i]*ref);
        		
        	}
            
    	
        }
        
        System.out.println();
        

	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 1.0;
        public int log10Steps = 6; // steps = (long)Math.pow(10, 6);
        
    }
}
