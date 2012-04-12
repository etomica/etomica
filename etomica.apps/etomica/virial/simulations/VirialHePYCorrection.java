package etomica.virial.simulations;

import etomica.api.IVectorMutable;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P2HeSimplified;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterDifference;
import etomica.virial.ClusterSum;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MayerXSpherical;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Computes corrections to Percus-Yevick approximations of B4 and B5 for a helium pair potential.
 * 
 * Select the compressibility or virial route via boolean variable compressibility.
 * 
 * Use semiclassical boolean to change pair potential to the quadratic Feynman-Hibbs potential.
 * 
 * Use calcApprox boolean to calculate quantities with the approximate potential.
 * Use subtractApprox boolean to compute the correction to that approximate potential.
 * 
 * If determining which option is most efficient via short calculations to estimate standard error, 
 * maintain a 50-50 split of steps between reference and target during data collection.
 * 
 * @author kate
 * @author Andrew Schultz
 */
public class VirialHePYCorrection {


    public static void main(String[] args) {

        VirialHePYCParam params = new VirialHePYCParam();
        boolean isCommandline = args.length > 0;
        ParseArgs.doParseArgs(params, args);
        
    	final int nPoints = params.nPoints;
    	final boolean compressibility = params.compressibility;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef < 0) {
            sigmaHSRef = 3 + 20/(10+temperatureK);
        }
        final boolean semiClassical = params.semiClassical;
        final double refFrac = params.refFrac;
        final boolean subtractApprox = params.subtractApprox;
        final boolean calcApprox = !subtractApprox && params.calcApprox;

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);

        System.out.println("Overlap sampling for He pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        if (semiClassical) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB[nPoints]);
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (subtractApprox) {
            System.out.println("computing difference from approximate He");
        }
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        MayerGeneralSpherical fTarget;
        MayerGeneralSpherical fTargetApprox;
        MayerXSpherical xTarget;
        MayerXSpherical xTargetApprox;
        if (semiClassical) {
            P2HeSimplified p2cApprox = new P2HeSimplified(space);
            Potential2Spherical p2Approx = p2cApprox.makeQFH(temperature);
            
            P2HePCKLJS p2c = new P2HePCKLJS(space);
            Potential2Spherical p2 = p2c.makeQFH(temperature);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);

            if (!compressibility) {
	    	    // need to add du to the semiclassical potential to do this.
	    	    throw new RuntimeException("can't do virial route with semiclassical potential");
	    	}
	    	xTarget = null;
            xTargetApprox = null;
        } else {
            P2HeSimplified p2Approx = new P2HeSimplified(space);
            
            P2HePCKLJS p2 = new P2HePCKLJS(space);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);

            xTarget = new MayerXSpherical(p2);
            xTargetApprox = new MayerXSpherical(p2Approx);
        }

        
        ClusterSum fullTargetCluster; 
        
        if (nPoints == 4) {
        	
        	if (compressibility) {
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0}}, {{0,2},{1,3}}};
    	        
    	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B4-B4PY(c)");
        		
    	        double[] weights = {-1.0/8.0};
    	     
    		    fullTargetCluster =  new ClusterSum(new ClusterBonds[] {cluster},weights,new MayerFunction[]{fTarget});
    		    
        	} else {
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B4-B4PY(v)");
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2},{1,3}}, {}};
    	        
    	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		int[][][] bondList2 = {{{0,2},{0,3},{1,2},{1,3}}, {{0,1}}};
        		
        		ClusterBonds cluster2 = new ClusterBonds(4, bondList2, false);
        		
        		double[] weights = {-1.0/8.0,1.0/12.0};
       	     
    		    fullTargetCluster =  new ClusterSum(new ClusterBonds[] {cluster, cluster2},weights,new MayerFunction[]{fTarget, xTarget});
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
	        	
	        	fullTargetCluster =  new ClusterSum(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4},weights,new MayerFunction[]{fTarget});
         
        	} else {
        		
        		System.out.println("Target diagram: Correction to Percus-Yevick approximation of virial coefficient: B5-B5PY(v)");
        		
        		int[][][] bondList1 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,3},{1,3},{2,4}}, {}, {{0,2},{1,4}} };
    	        
    	        ClusterBonds c1 = new ClusterBonds(5, bondList1, false);
    	        
    	        int[][][] bondList2 = { {{1,2},{2,3},{3,4},{4,0},{0,2},{1,4}}, {}, {{0,1},{0,3},{1,3},{2,4}} };
    	        
    	        ClusterBonds c2 = new ClusterBonds(5, bondList2, false);
    	        
    	        int[][][] bondList3 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{0,2},{0,3},{1,3},{1,4},{2,4}}, {} };
    	        
    	        ClusterBonds c3 = new ClusterBonds(5, bondList3, false);
    	        
    	        int[][][] bondList4 = { {{0,1},{1,2},{2,3},{3,4},{4,0}}, {{2,4}} };
    	        
    	        ClusterBonds c4 = new ClusterBonds(5, bondList4, false);
    	        
    	        int[][][] bondList5 = { {{0,1},{1,2},{2,3},{3,4},{4,0},{1,3}}, {{0,3}} };
    	        
    	        ClusterBonds c5 = new ClusterBonds(5, bondList5, false);
    	        
    	    	double[] weights = {0.5,-1.0/3.0,-0.2, 1.0/6.0, 1.0/3.0};
          	     
    		    fullTargetCluster =  new ClusterSum(new ClusterBonds[] {c1, c2, c3, c4, c5},weights,new MayerFunction[]{fTarget, xTarget});

        	}
        	
        } else {
        	
        	throw new IllegalArgumentException("Cannot yet compute correction to Percus-Yevick approximation for that order of virial coefficient");
        	
        }

        ClusterAbstract targetCluster = null;
        if (subtractApprox) {
            final ClusterSum[] targetSubtract = new ClusterSum[1];
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            if (compressibility) {
                targetSubtract[0] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetApprox});
            }
            else {
                targetSubtract[0] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetApprox, xTargetApprox});
            }
            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
        }
        else {
            targetCluster = fullTargetCluster;
        }

        
        targetCluster.setTemperature(temperature);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("He")), temperature, refCluster, targetCluster);
        sim.integratorOS.setAgressiveAdjustStepFraction(true);
        
        long t1 = System.currentTimeMillis();
        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        double r = 4;
        for (int i=1; i<nPoints; i++) {
            IVectorMutable v = sim.box[1].getLeafList().getAtom(i).getPosition();
            v.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
        if (sim.box[1].getSampleCluster().value(sim.box[1]) == 0) {
            throw new RuntimeException("oops");
        }
        
        sim.integratorOS.setNumSubSteps(1000);
        
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        steps /= 1000;
        sim.setAccumulatorBlockSize(steps);
        
        System.out.println();
        String refFileName = null;
        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints+"_2b_"+tempString;
            refFileName += semiClassical ? "_sc" : "_c";
            if (calcApprox) {
                refFileName += "a";
            }
            else if (subtractApprox) {
                refFileName += "sa";
            }
            refFileName += "PY";
            if (compressibility) {
                refFileName += "C";
            }
            else {
                refFileName += "V";
            }
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().actionPerformed();

        System.out.println();
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println();
        
        sim.printResults(HSB[nPoints]);

        DataGroup allYourBase = (DataGroup)sim.accumulators[1].getData();
        IData averageData = allYourBase.getData(sim.accumulators[1].AVERAGE.index);
        IData errorData = allYourBase.getData(sim.accumulators[1].ERROR.index);
        IData covarianceData = allYourBase.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        int n = 0;
        double correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("diagram "+nPoints+"bc average: %20.15e error: %9.4e ocor: %6.4f\n",
                averageData.getValue(0), errorData.getValue(0), correlationCoef));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHePYCParam extends ParameterBase {
        public int nPoints = 4;
        public boolean compressibility = true;
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = -1;
        public boolean semiClassical = true;
        public boolean calcApprox = false;
        public boolean subtractApprox = true;
    }
}
