/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * Computes corrections to hypernetted-chain approximations of B4 and B5 for a helium pair potential.
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
public class VirialHePCKLJSCHNCCorrection {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();

        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; boolean compressibility; int stepsPerBlock; long eqSteps;
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
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef +" Angstroms");
        
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
        
        MayerESpherical eTarget;
        if (QFH) {
        	P2HePCKLJS p20 = new P2HePCKLJS(space);
        	P2EffectiveFeynmanHibbs p2 = new P2EffectiveFeynmanHibbs(Space3D.getInstance(),p20);
			p2.setTemperature(temperature);		
			p2.setMass(4.002602);
	    	eTarget = new MayerESpherical(p2);
        } else {
        	P2HePCKLJS p2 = new P2HePCKLJS(space);
        	eTarget = new MayerESpherical(p2);
        }
            
        ClusterAbstract targetCluster;
        
        if (nPoints == 4) {
        	
        	if (compressibility) {
        	
        		System.out.println("Target diagram: Correction to HNC(c) approximation of virial coefficient: B4-B4HNC(c)");
        
	        	int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2}}, {{1,3}} };
	        
	        	ClusterBonds cluster = new ClusterBonds(4, bondList, false);
	
	        	double[] weights = {-1.0/8.0};
		    
		    	targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eTarget});
		    
        	} else {
        		
        		System.out.println("Target diagram: Correction to HNC(v) approximation of virial coefficient: B4-B4HNC(v)");
        		
        		int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2},{1,3}}, {}};

        		ClusterBonds cluster = new ClusterBonds(4, bondList, false);
        		
        		double[] weights = {-1.0/8.0};
		    
        		targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eTarget});
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
	        	
	        	targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4, cluster5, cluster6},weights,new MayerFunction[]{eTarget});
        	
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
	        	
	        	targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4},weights,new MayerFunction[]{eTarget});
        	}
        	
        } else {
        	
        	throw new IllegalArgumentException("Cannot yet compute correction to CHNC(c) approximation for that order of virial coefficient");
        	
        }

        
        
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" IntegratorOverlap steps of 1000 steps)");
		
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactorySpheres(), temperature,refCluster,targetCluster);
        
        
        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box[1]);
        
        sim.integratorOS.setNumSubSteps(1000);
        
        // Keep 50-50 split between reference and target
        sim.integratorOS.setAdjustStepFreq(adjustStepFrequency);
        sim.setAccumulatorBlockSize(stepsPerBlock);
        System.out.println(stepsPerBlock+" steps per block");
        
        System.out.println();
        
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        sim.initRefPref(refFileName, steps/100);

        
        sim.equilibrate(refFileName, eqSteps); // 5000 IntegratorOverlap steps = 5e6 steps
        System.out.println((eqSteps*1000) + " equilibration steps (" + eqSteps + " Integrator Overlap Steps)"); 
        
        System.out.println("equilibration finished");

        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().actionPerformed();
        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());

        System.out.println();
        
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);

        System.out.println();
        System.out.println("cm"+((nPoints-1)*3)+"/mol"+(nPoints-1)+": ");
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1)+", error: "+ratioAndError[1]*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1));
        System.out.println();
        System.out.println("Simulation units:");
        
        System.out.println();
        
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("reference ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].RATIO_ERROR.index)).getData()[1]);
        System.out.println("reference average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[0]);
        System.out.println("reference overlap average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);
        
        
        System.out.println("reference autocorrelation function: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].BLOCK_CORRELATION.index)).getData()[0]);
        
        System.out.println("reference overlap autocorrelation function: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].BLOCK_CORRELATION.index)).getData()[1]);
        
        System.out.println();
        
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].RATIO_ERROR.index)).getData()[1]);
        System.out.println("target average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].ERROR.index)).getData()[0]);
        System.out.println("target overlap average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[1].ERROR.index)).getData()[1]);

        System.out.println("target autocorrelation function: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].BLOCK_CORRELATION.index)).getData()[0]);
        
        System.out.println("target overlap autocorrelation function: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].BLOCK_CORRELATION.index)).getData()[1]);
    
    }

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
    	
    	// number of molecules in simulation (e.g., 2 for B2 calculation)

 
        public int numMolecules = 4;
        public boolean compressibility = false;
        public double temperature = 5;
        public long numSteps = 1000;
        public int stepsPerBlock = 1000;
        public long eqSteps = 1000;
        public double sigmaHSRef = 3;
        public boolean writeRefPref = false;
        public boolean QFH = true;
        public boolean adjustStepFrequency = false;
    }
}
