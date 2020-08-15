/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;




import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationClusterMove;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerFunction;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * 
 * An implementation of MSMC with direct sampling for computation of 
 * virial coefficients of hard spheres, in the spirit of Ree-Hoover Monte Carlo 
 * (Ree & Hoover, 1964, Fifth and Sixth Virial Coefficients for Hard Spheres and 
 * Hard Disks).  The reference value is not a virial coefficient.
 * 
 * Kate Shaul (2009)
 */

public class ReeHooverMC {


    public static void main(String[] args) {

    	VirialParam params = new VirialParam();
 
        final int n; int log10Steps; double sigmaHSRef;
        
        if (args.length == 0) {
        	
        	n = params.n;
        	sigmaHSRef = params.sigmaHSRef;
            log10Steps = params.log10Steps;
            
        } else if (args.length == 2) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	n = Integer.parseInt(args[0]);
        	sigmaHSRef = Double.parseDouble(args[1]);
            log10Steps = Integer.parseInt(args[2]);

        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments.");
        }

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
        System.out.println("B8HS: "+HSB[8]+" = 0.004164 B2HS^7");

        
        long steps = (long)Math.pow(10, log10Steps);

        double b = 2*Math.PI/3;
        
        double ref = Math.pow(-2*b, n-1);

        
        
        System.out.println("HS B" + n);
       
        System.out.println();
        System.out.println("Direct sampling with *presence* of linear chain as reference, in accordance with Ree-Hoover MC method.");
        System.out.println("Reference value for " + n + " = " + ref);
        System.out.println();
        System.out.println("10^" + log10Steps + " steps");
        System.out.println();
		
        Space space = Space3D.getInstance();
        
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        
        ClusterAbstract[] targetCluster = new ClusterAbstract[1];
        targetCluster[0]= Standard.virialCluster(n, fRef, n>3, eRef, true);
        targetCluster[0].setTemperature(1.0);

        ClusterAbstract refCluster;
        
        if (n == 4) {
        	
        	int[][][] bondList = {{{0,1},{1,2},{2,3}}, {}};
 	        
 	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
 	
 	        double[] weights = {1.0};
 	  
 		    
 		    refCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eRef});
        }
        
        else if (n == 5) {
        	 int[][][] bondList = {{{0,1},{1,2},{2,3}, {3,4} }, {}};
 	        
 	        ClusterBonds cluster = new ClusterBonds(5, bondList, false);
 	        
 	        double[] weights = {1.0};
 	  
 		    
 		    refCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eRef});
        } else {
        	
        	throw new IllegalArgumentException("This method only work for n = 4 and 5.");
        
	       
        }
        
        refCluster.setTemperature(1.0);
        
        //////////////////////////////////////////////////////////////////////////////////////////////
        // The sampleCluster is the system used to sample configuration space.  This should be the full LJ potential, as its important configuration space 
        // encompasses that of LJTS.
        //////////////////////////////////////////////////////////////////////////////////////////////
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        
        SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
		double temp = 1.0;
        
        final SimulationVirial sim = new SimulationVirial(space, speciesFactory, temp, sampleCluster,refCluster,targetCluster);
        
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box);
        
        
        
        sim.equilibrate(steps/40);
        
        System.out.println("equilibration finished");
        
        System.out.println();

        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
sim.getController2().runActivityBlocking(new etomica.action.activity.ActivityIntegrate2(sim.integrator), steps);
        
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        
        System.out.println();
        System.out.println("Linear chain: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
                           +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
                           +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);
        
        double ratio;
        double error;
        
        System.out.println();
        System.out.println("HS");
        System.out.println();
        
    	
    	ratio = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];
        error = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];
        
    	System.out.println("HS average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1]
                               +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[1]
                               +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1]);
       
        
    	
    	System.out.println();
        

        System.out.println(" ratio average: "+ratio+", error: "+error);

        
        System.out.println();
        
        		
        System.out.println(" Bn: "+ratio*ref+", error: "+ error*ref);

        
        System.out.println();
        

	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int n = 5;
        public double sigmaHSRef = 1;
        public int log10Steps = 8; // steps = (long)Math.pow(10, 6);

    }
}
