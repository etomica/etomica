/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * LJ simulation using Mayer sampling to evaluate corrections to the HNC(c) approximations to B4 and B5
 */
public class VirialLJCHNCCorrection {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();

        double temperature; final int nPoints; double sigmaHSRef;
        long steps; boolean compressibility; int stepsPerBlock; long eqSteps;
        if (args.length == 0) {
        	
        	nPoints = params.numMolecules;
        	compressibility = params.compressibility;
            temperature = params.temperature;
            steps = params.numSteps;
            stepsPerBlock = params.stepsPerBlock;
            eqSteps = params.eqSteps;
            sigmaHSRef = params.sigmaHSRef;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 7) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	compressibility = Boolean.parseBoolean(args[1]);
        	temperature = Double.parseDouble(args[2]);
            steps = Integer.parseInt(args[3]);
            stepsPerBlock = Integer.parseInt(args[4]);
            eqSteps = Integer.parseInt(args[5]);
            sigmaHSRef = Double.parseDouble(args[6]);
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

        
        System.out.println("Overlap sampling for Lennard-Jones system at reduced temperature " + temperature);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef);
        
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
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        
        ClusterAbstract targetCluster;
        
        
        
        if (nPoints == 4) {
        	
        	if (compressibility) {
        	
        		System.out.println("Target diagram: Correction to HNC(c) approximation of virial coefficient: B4-B4HNC(c)");
        
	        	int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2}}, {{1,3}} };
	        
	        	ClusterBonds cluster = new ClusterBonds(4, bondList, false);
	
	        	double[] weights = {-1.0/8.0};
		    
		    	targetCluster =  new ClusterSum(new ClusterBonds[] {cluster},weights,new MayerFunction[]{fTarget});
		    
        	} else {
        		
        		System.out.println("Target diagram: Correction to HNC(v) approximation of virial coefficient: B4-B4HNC(v)");
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2},{1,3}}, {}};
        		
        		ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		double[] weights = {-1.0/8.0};
		    
        		targetCluster =  new ClusterSum(new ClusterBonds[] {cluster},weights,new MayerFunction[]{fTarget});
        	}
		    
        } else if (nPoints == 5) {
        	
        	
        	if (compressibility) {
        	
	        	System.out.println("Target diagram: Correction to HNC(c) approximation of virial coefficient: B5-B5HNC(c)");
	        	
	        	int[][][] bondList1 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{2,4}}, {{0,2},{0,3},{1,3},{1,4}} };
	        	
	        	int[][][] bondList2 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{2,4}}, {{0,3},{1,3},{1,4}} };
	        	
	        	int[][][] bondList3 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{1,4}}, {{0,3},{1,3},{2,4}} };
	        	
	        	int[][][] bondList4 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,3},{1,3},{2,4}}, {{0,2},{1,4}} };
	        	
	        	int[][][] bondList5 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{0,3},{1,3},{1,4}}, {{2,4}} };
	        	
	        	int[][][] bondList6 = { {{1,2},{2,3},{3,4},{4,0},{0,2},{1,4},}, {{0,1},{0,3},{1,3},{2,4}} };
	        	
	        	ClusterBonds cluster1 = new ClusterBonds(5, bondList1, false);
	        	ClusterBonds cluster2 = new ClusterBonds(5, bondList2, false);
	        	ClusterBonds cluster3 = new ClusterBonds(5, bondList3, false);
	        	ClusterBonds cluster4 = new ClusterBonds(5, bondList4, false);
	        	ClusterBonds cluster5 = new ClusterBonds(5, bondList5, false);
	        	ClusterBonds cluster6 = new ClusterBonds(5, bondList6, false);
	        	
	        	double[] weights = {-0.4, 0.2, 0.1, 0.6, -13.0/30.0, -0.1};
	        	
	        	targetCluster =  new ClusterSum(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4, cluster5, cluster6},weights,new MayerFunction[]{fTarget});
        	
        	} else {
	        
        		System.out.println("Target diagram: Correction to HNC(v) approximation of virial coefficient: B5-B5HNC(v)");
	        	
	        	int[][][] bondList1 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{1,4}}, {{0,3},{1,3},{2,4}} };
	        	
	        	int[][][] bondList2 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,3},{1,3},{2,4}}, {{0,2},{1,4}} };
	        	
	        	int[][][] bondList3 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{0,3},{1,3},{1,4}}, {{2,4}} };
	        	
	        	int[][][] bondList4 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{0,3},{1,3},{1,4},{2,4}}, {} };
	        	
	        	ClusterBonds cluster1 = new ClusterBonds(5, bondList1, false);
	        	ClusterBonds cluster2 = new ClusterBonds(5, bondList2, false);
	        	ClusterBonds cluster3 = new ClusterBonds(5, bondList3, false);
	        	ClusterBonds cluster4 = new ClusterBonds(5, bondList4, false);
	        	
	        	double[] weights = {-1.0, 1.5, -1.0/3.0, -0.2};
	        	
	        	targetCluster =  new ClusterSum(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4},weights,new MayerFunction[]{fTarget});
        	}
        	
        } else {
        	
        	throw new IllegalArgumentException("Cannot yet compute correction to CHNC(c) approximation for that order of virial coefficient");

        }


        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints > 3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps * 1000) + " steps (" + steps + " IntegratorOverlap steps of 1000 steps)");

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A"))), nPoints, temperature, refCluster, targetCluster);


        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box[1]);

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAdjustStepFraction(false);
        sim.setAccumulatorBlockSize(stepsPerBlock);
        System.out.println(stepsPerBlock + " steps per block");
        
        System.out.println();
        
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        /*final long steps2=steps;
        IAction progressReport = new IAction() {
            public void actionPerformed() {
                
            	 double stepSizeDivider = ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).getAdjustStepSize();
            	 
            	 double stepSize = sim.mcMoveTranslate[1].getStepSize();
            	 
            	
            	 double rate = ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).acceptanceRatio();
            	 double stepSizeChange = stepSize*(stepSizeDivider-1);
            	 
            	 System.out.println(stepSize + "  " + stepSizeChange + "  " + rate);
            	 
            	 if ( (Math.abs(rate-0.5) < 0.001) && (sim.integratorOS.getStepCount() > (steps2/40)) ){
            		 System.out.println((sim.integratorOS.getStepCount()*1000) + " equilibration steps"); 
            		 sim.ai.setMaxSteps(0);
            		 
            	 }
            	 
            }
        };   
        
        IntegratorListenerAction iLA = new IntegratorListenerAction(progressReport, (int)(100));
        sim.integratorOS.getEventManager().addListener(iLA);
		*/
        
        sim.equilibrate(refFileName, eqSteps); // 5000 IntegratorOverlap steps = 5e6 steps
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, steps);
System.out.println((eqSteps*1000) + " equilibration steps (" + eqSteps + " Integrator Overlap Steps)");

        //sim.integratorOS.getEventManager().removeListener(iLA);

        System.out.println("equilibration finished");

    /*    IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double ratio = sim.dsvo.getDataAsScalar();
                double error = sim.dsvo.getError();
                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
            }
        };
        sim.integratorOS.addIntervalAction(progressReport);
        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));*/

        for (int i = 0; i < 2; i++) {
            System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(ai);
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
        public boolean compressibility = false;
        public double temperature = 0.6;
        public long numSteps = 10000;
        public int stepsPerBlock = 1000;
        public long eqSteps = 1000;
        public double sigmaHSRef = 1.5;
        public boolean writeRefPref = false;
    }
}
