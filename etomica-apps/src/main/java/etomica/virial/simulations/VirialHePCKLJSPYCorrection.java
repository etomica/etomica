/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * Computes corrections to Percus-Yevick approximations of B4 and B5 for a helium pair potential.
 * 
 * Select the compressibility or virial route via boolean variable compressibility.
 * 
 * Use QFH boolean to change pair potential to the quadratic Feynman-Hibbs potential.
 * 
 * If determining which option is most efficient via short calculations to estimate standard error, 
 * maintain a 50-50 split of steps between reference and target during data collection with 
 * with adjustStepFreqFrequency = false;
 * 
 * @author kate
 *
 */
public class VirialHePCKLJSPYCorrection {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();

        boolean compressibility;
        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; int stepsPerBlock; long eqSteps;
        boolean QFH; boolean adjustStepFrequency;
        if (args.length == 0) {
        	
        	nPoints = params.numMolecules;
        	compressibility = params.compressibility;
            temperatureK = params.temperature;
            steps = params.numSteps;
            stepsPerBlock = params.stepsPerBlock;
            eqSteps = params.eqSteps;
            sigmaHSRef = params.sigmaHSRef;
            QFH = params.QFH;
            adjustStepFrequency = params.adjustStepFrequency;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 9) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	compressibility = Boolean.parseBoolean(args[1]);
        	temperatureK = Double.parseDouble(args[2]);
            steps = Integer.parseInt(args[3]);
            stepsPerBlock = Integer.parseInt(args[4]);
            eqSteps = Integer.parseInt(args[5]);
            sigmaHSRef = Double.parseDouble(args[6]);
            QFH = Boolean.parseBoolean(args[7]);
            adjustStepFrequency = Boolean.parseBoolean(args[8]);
            params.writeRefPref = true;
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        } 

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);

        System.out.println("Overlap sampling for He pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        if (QFH) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B2HS: "+HSB[2]);
        System.out.println("  B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("  B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("  B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("  B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("  B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
        System.out.println("  B8HS: "+HSB[8]+" = 0.004164 B2HS^7");
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        
        MayerGeneralSpherical fTarget; MayerESpherical eTarget; MayerXSpherical xTarget;
        if (QFH) {
        	P2HePCKLJS p20 = new P2HePCKLJS(space);
        	P2EffectiveFeynmanHibbs p2 = new P2EffectiveFeynmanHibbs(Space3D.getInstance(),p20);
			p2.setTemperature(temperature);		
			p2.setMass(4.002602);
			fTarget = new MayerGeneralSpherical(p2);
	    	eTarget = new MayerESpherical(p2);
	    	xTarget = new MayerXSpherical(p20);
        } else {
        	P2HePCKLJS p2 = new P2HePCKLJS(space);
        	fTarget = new MayerGeneralSpherical(p2);
        	eTarget = new MayerESpherical(p2);
        	xTarget = new MayerXSpherical(p2);
        }

        
        ClusterAbstract targetCluster; 
        
        if (nPoints == 4) {
        	
        	if (compressibility) {
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0}}, {{0,2},{1,3}}};
    	        
    	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B4-B4PY(c)");
        		
    	        double[] weights = {-1.0/8.0};
    	     
    		    targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eTarget});
    		    
        	} else {
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B4-B4PY(v)");
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2},{1,3}}, {}, {}};
    	        
    	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		int[][][] bondList2 = {{{0,2},{0,3},{1,2},{1,3}}, {}, {{0,1}}};
        		
        		ClusterBonds cluster2 = new ClusterBonds(4, bondList2, false);
        		
        		double[] weights = {-1.0/8.0,1.0/12.0};
       	     
    		    targetCluster =  new ClusterSum(new ClusterBonds[] {cluster, cluster2},weights,new MayerFunction[]{fTarget, eTarget, xTarget});
        	}
		    
        } else if (nPoints == 5) { 
        	
        	if (compressibility) {
        		
        		// B5-B5PY = -1/5*S1 + 1/2*S2 - 1/3*S3 + S4 
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B5-B5PY(c)");
        	
	        	int[][][] bondList1 = { {{0,1},{1,2},{2,3},{3,4},{4,0}}, {{0,2},{0,3},{1,3},{1,4},{2,4}} };
		        
		        ClusterBonds cluster1 = new ClusterBonds(5, bondList1, false);
		        
		        int[][][] bondList2 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,3},{1,3},{2,4}}, {{0,2},{1,4}} };
		        
		        ClusterBonds cluster2 = new ClusterBonds(5, bondList2, false);
		        
		        int[][][] bondList3 = { {{1,2},{2,3},{3,4},{4,0},{0,2},{1,4}}, {{0,1},{0,3},{1,3},{2,4}} };
		        
		        ClusterBonds cluster3 = new ClusterBonds(5, bondList3, false);
		        
		        int[][][] bondList4 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{1,3},{2,4}}, {{0,2},{1,4}} };
		        
		        ClusterBonds cluster4 = new ClusterBonds(5, bondList4, false);
	        	
		        double[] weights = {-0.2, 0.5, -1.0/3.0, 1.0};
	        	
	        	targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4},weights,new MayerFunction[]{eTarget});
         
        	} else {
        		
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B5-B5PY(v)");
        		
        		int[][][] bondList1 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,3},{1,3},{2,4}}, {{0,2},{1,4}} ,{}};
    	        
    	        ClusterBonds c1 = new ClusterBonds(5, bondList1, false);
    	        
    	        int[][][] bondList2 = { {{1,2},{2,3},{3,4},{4,0},{0,2},{1,4}}, {{0,1},{0,3},{1,3},{2,4}}, {} };
    	        
    	        ClusterBonds c2 = new ClusterBonds(5, bondList2, false);
    	        
    	        int[][][] bondList3 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{0,3},{1,3},{1,4},{2,4}}, {}, {} };
    	        
    	        ClusterBonds c3 = new ClusterBonds(5, bondList3, false);
    	        
    	        int[][][] bondList4 = { {{0,1},{1,2},{2,3},{3,4},{4,0}}, {}, {{2,4}} };
    	        
    	        ClusterBonds c4 = new ClusterBonds(5, bondList4, false);
    	        
    	        int[][][] bondList5 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{1,3}}, {}, {{0,3}} };
    	        
    	        ClusterBonds c5 = new ClusterBonds(5, bondList5, false);
    	        
    	    	double[] weights = {0.5,-1.0/3.0,-0.2, 1.0/6.0, 1.0/3.0};
          	     
    		    targetCluster =  new ClusterSum(new ClusterBonds[] {c1, c2, c3, c4, c5},weights,new MayerFunction[]{fTarget, eTarget, xTarget});

        	}
        	
        } else {
        	
        	throw new IllegalArgumentException("Cannot yet compute correction to Percus-Yevick approximation for that order of virial coefficient");
        	
        }

        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" IntegratorOverlap steps of 1000 steps)");

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")), temperature, refCluster, targetCluster);
        
        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box[1]);
        
        sim.integratorOS.setNumSubSteps(1000);
        
        // Keep 50-50 split between reference and target
        sim.integratorOS.setAdjustStepFraction(adjustStepFrequency);
       
        sim.setAccumulatorBlockSize(stepsPerBlock);
        System.out.println(stepsPerBlock+" steps per block");
        
        System.out.println();
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;

        sim.initRefPref(refFileName, steps/100);
        sim.equilibrate(refFileName, eqSteps); // 5000 IntegratorOverlap steps = 5e6 steps
        System.out.println((eqSteps*1000) + " equilibration steps (" + eqSteps + " Integrator Overlap Steps)");

        System.out.println("equilibration finished");


        for (int i = 0; i < 2; i++) {
            System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), steps);

        System.out.println();
        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

        System.out.println();

        sim.printResults(HSB[nPoints]);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
    	
    	// number of molecules in simulation (e.g., 2 for B2 calculation)

 
        public int numMolecules = 4;
        public boolean compressibility = true;
        public double temperature = 500;
        public long numSteps = 1000;
        public int stepsPerBlock = 1000;
        public long eqSteps=1000;
        public double sigmaHSRef = 3;
        public boolean writeRefPref = false;
        public boolean QFH = true;
        public boolean adjustStepFrequency = false;
        
    }
}
