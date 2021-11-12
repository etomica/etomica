/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.squarewell;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerFunction;
import etomica.virial.MeterVirialSWWE;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterAtomHSChain;
import etomica.virial.mcmove.MCMoveClusterAtomHSRing;
import etomica.virial.mcmove.MCMoveClusterAtomHSTree;
import etomica.virial.simulations.SimulationVirial;
import etomica.virial.simulations.hardsphere.VirialHS.VirialHSParam;
import etomica.virial.wheatley.ClusterWheatleyExtendSW;


public class VirialSWWE {
	
	public static void main(String [] args){
		
		VirialSWWEParam params = new VirialSWWEParam();
		if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }else {
            params.nPoints = 4;
            params.numSteps = 100000L;
            params.ref = VirialHSParam.RING_CHAIN_TREES;
            params.chainFrac = 0.1;
            params.treeFrac = 0.1;
            params.ringFrac = (1- params.chainFrac - params.treeFrac);
        }
		
		final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;
        final double sigmaHS = 2.0;
        final double sigmaSW = 1.0;
        final double lambda = 2.0;
        final double chainFrac = params.chainFrac;
        final double treeFrac = params.treeFrac;
        final double ringFrac = (1 - chainFrac - treeFrac);
        
        double vhs = (4.0/3.0)*Math.PI*sigmaHS*sigmaHS*sigmaHS;
        
        Space space = Space3D.getInstance();
        
        MayerEHardSphere fTargete2 = new MayerEHardSphere(1.0);
        MayerFunction fTargetf1 = new MayerFunction() {
            final double sigma2 = 1.0;
            final double well2 = lambda*lambda;
            
            public void setBox(Box box) {}
            
            public IPotential getPotential() {return null;}
            
            public double f(IMoleculeList pair, double r2, double beta) {
                if (r2 < sigma2 || r2 > well2) return 0;
                return 1;
            }
        };
        MayerFunction fRefPos = new MayerFunction() {

            public void setBox(Box box) {}
            public IPotential getPotential() {return null;}

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHS*sigmaHS ? 1 : 0;
            }
        };
        
        ClusterWheatleyExtendSW targetCluster = new ClusterWheatleyExtendSW(nPoints, fTargetf1, fTargete2);
        
        targetCluster.setTemperature(1.0);

        ClusterAbstract refCluster = null;
        long numDiagrams = 0;
        
        if (ref == VirialSWWEParam.RING_CHAIN_TREES){
        	/* This is from VirialBinMultiThreaded.java, it is efficient.*/
        	ClusterChainHS cr = new ClusterChainHS(nPoints, fRefPos, true);//ring
            long numRingDiagrams = cr.numDiagrams();
            
            double ringIntegral = numRingDiagrams*Standard.ringHS(nPoints)*Math.pow(sigmaHS, 3*(nPoints-1));//ring integral: Coeffi*(sigmaHS^3)^(N-1)
//            System.out.println("---ring integral: "+ringIntegral);
            double chainIntegral = (SpecialFunctions.factorial(nPoints)/2)*Math.pow(vhs, nPoints-1);//N!/2*(4/3*pi*sigma^3)^(N-1)
//            System.out.println("---chain integral: "+chainIntegral);
            
            ClusterChainHS crc = new ClusterChainHS(nPoints, fRefPos, chainFrac/chainIntegral, ringFrac/ringIntegral);//chain and ring
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);//tree
            
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{crc, ct});
            
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }
            double treeIntegral = numTreeDiagrams*Math.pow(vhs, nPoints-1);//N^(N-2)*(4/3*pi*sigma^3)^(N-1)
//            System.out.println("---tree integral: "+treeIntegral);
            // weighting for chain and ring are handled internally
            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{1,treeFrac/treeIntegral});
//          System.out.println("Inside if1---> chainFrac = " + chainFrac + ", treeFrac = " + treeFrac + ", ringFrac = " + ringFrac);

            
            
            /* This is from VirialHS.java, it is not efficient.*/
//        	ClusterChainHS cr = new ClusterChainHS(nPoints, fRefPos, true);//ring
//        	long numRingDiagrams = cr.numDiagrams();
//        	ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);//chain
//        	ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);//tree
//        	
//        	refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cr, cc, ct});
//        	
//        	long numTreeDiagrams = 1;
//            for (int i=0; i<nPoints-2; i++) {
//                numTreeDiagrams *= nPoints;
//            }
            
//            final double dr = 0.00001;
//            CalcFFT myFFT = new CalcFFT(new IFunction() {
//                public double f(double x) {
//                    if (Math.abs(x-1) < 0.1*dr) {
//                        return 0.5;
//                    }
//                    return x<1 ? 1 : 0;
//                }
//            }, dr, 20);
//            
//            List<Object> strands = new ArrayList<Object>();
//            strands.add(2);
//            List<Integer> list1 = new ArrayList<Integer>();
//            for (int i=1; i<nPoints; i++) {
//                list1.add(2);
//            }
//            strands.add(list1);
//            List<Object> oneMore = new ArrayList<Object>();
//            oneMore.add(strands);
            
//            double ringIntegral = numRingDiagrams*Standard.ringHS(nPoints);
//            System.out.println("---ring integral: "+ringIntegral);
//            double chainIntegral = (SpecialFunctions.factorial(nPoints)/2)*Math.pow(vhs, nPoints-1);//N!/2*(4/3*pi*sigma^3)^(N-1)
//            System.out.println("---chain integral: "+chainIntegral);
//            double treeIntegral = numTreeDiagrams*Math.pow(vhs, nPoints-1);//N^(N-2)*(4/3*pi*sigma^3)^(N-1)
//            System.out.println("---tree integral: "+treeIntegral);
            
//            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{ringFrac/ringIntegral,chainFrac/chainIntegral,treeFrac/treeIntegral});
//            System.out.println("Inside if1---> chainFrac = " + chainFrac + ", treeFrac = " + treeFrac + ", ringFrac = " + ringFrac);	
        }
        
        final double Coeff = 1;
        refCluster.setTemperature(1.0);
        
        ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};

        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")));
        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{species}, new int[]{nPoints}, 1.0, ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, targetDiagrams);
        sim.setDoFasterer(true);
        sim.setDoWiggle(false);
        sim.init();
        MeterVirialSWWE meter = new MeterVirialSWWE(targetCluster);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        
        AccumulatorAverageCovariance accumulator = new AccumulatorAverageCovariance(1, true);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(100000000);
        
        sim.integratorFasterer.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        if (ref == VirialSWWEParam.RING_CHAIN_TREES){
        	MCMoveClusterAtomHSRing mcMoveHSR = new MCMoveClusterAtomHSRing(sim.getRandom(), sim.box, sigmaHS);
            sim.integratorFasterer.getMoveManager().addMCMove(mcMoveHSR);
            sim.integratorFasterer.getMoveManager().setFrequency(mcMoveHSR, ringFrac);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), sim.box, sigmaHS);
            sim.integratorFasterer.getMoveManager().addMCMove(mcMoveHSC);
            sim.integratorFasterer.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), sim.box, sigmaHS);
            sim.integratorFasterer.getMoveManager().addMCMove(mcMoveHST);
            sim.integratorFasterer.getMoveManager().setFrequency(mcMoveHST, treeFrac);
            
//            System.out.println("Inside if2---> chainFrac = " + chainFrac + ", treeFrac = " + treeFrac + ", ringFrac = " + ringFrac);
        }
        System.out.println("\n**********Calculate B" + nPoints +" of square well potential************\n");
        System.out.println("	Display input arguments:");
        System.out.println("	# of total steps = " + steps);
        System.out.println("	sigmaHSRef = " + sigmaHS );
        System.out.println("	sigmaSW = " + 1);
        System.out.println("	lamda = " + lambda);
        System.out.println("	Y_Value = " + 1);
        System.out.println("	chainFrac = " + chainFrac + ", treeFrac = " + treeFrac + ", ringFrac = " + ringFrac);

        System.out.println("\n	*****" + steps + " steps for chain, tree and ring generation ");


        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorFasterer, steps));
        long t2 = System.currentTimeMillis();
        
        
        DataGroup allYourBase = (DataGroup)accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);
        IData covarianceData = allYourBase.getData(accumulator.COVARIANCE.index);
        
        System.out.println();
        
        System.out.println("	coefficients = ");
        for(int i=0; i<averageData.getLength(); i++){
        	double avg = averageData.getValue(i);
        	double err = errorData.getValue(i);	
//        	System.out.print(String.format("	    %17.13e  %7.7e\n", avg*refIntegral, err*Math.abs(refIntegral)));
        	System.out.print(String.format("	    %17.13e  %7.7e\n", avg*Coeff, err*Math.abs(Coeff)));
        }	
        System.out.println("	expectation between coefficients = ");
        int k=0;
        for(int i=0; i<averageData.getLength(); i++){//Calculate the correlation from the covariance matrix.
			for(int j=i+1; j<averageData.getLength(); j++){
				double covariance = covarianceData.getValue(i*averageData.getLength()+j);
				double avgI = averageData.getValue(i);
				double avgJ = averageData.getValue(j);
//				System.out.println("covariance="+covariance+";");
//				System.out.println("avgI="+avgI+";");
//				System.out.println("avgJ="+avgJ+";");
				double expectation = covariance + avgI*avgJ;
				System.out.print(String.format("	   %3d   %25.15e between system[%d] and system[%d]\n", k, expectation, i, j));
				k++;
			}
		}
        
        System.out.println("	correlation between coefficients (sampling by block)= ");
        k=0;
        for(int i=0; i<averageData.getLength(); i++){//Calculate the correlation from the covariance matrix.
			for(int j=i+1; j<averageData.getLength(); j++){
				double covariance = covarianceData.getValue(i*averageData.getLength()+j);
				double varianceI = covarianceData.getValue(i*averageData.getLength()+i);
				double varianceJ = covarianceData.getValue(j*averageData.getLength()+j);
//				System.out.println("covariance="+covariance+";");
//				System.out.println("varianceI="+varianceI+";");
//				System.out.println("varianceJ="+varianceJ+";");
				double correlation;
				if (varianceI * varianceJ != 0){
					correlation = covariance / Math.sqrt(varianceI*varianceJ);
				}else{
					correlation = 0.0;
				}
				System.out.print(String.format("	   %3d   %25.15e between B%d[%d] and B%d[%d]\n", k, correlation, averageData.getLength(), i, averageData.getLength(), j));
				k++;
			}
		}
        
        System.out.print("\n	Time for " + steps + " steps = ");
        System.out.print(String.format(		"%.5f seconds = ", (t2-t1)/1000.0));
        System.out.print(String.format(		"%.5f minutes = ", (t2-t1)/1000.0/60.0));
        System.out.print(String.format(		"%.5f hours = ", (t2-t1)/1000.0/3600.0));
        System.out.println(String.format(	"%.5f days", (t2-t1)/1000.0/3600.0/24));
        
        System.out.print("\n	Speed = ");
        System.out.print(String.format(		"%.6f step/seconds = ", (double)steps/((t2-t1)/1000.0)));
        System.out.print(String.format(		"%.6f thousand_step/second = ", ((double)steps/((t2-t1)/1000.0))/1000.0));
        System.out.println(String.format(		"%.6f million_step/second", ((double)steps/((t2-t1)/1000.0))/1000000.0));
        System.out.println("\n**********End***********\n");
	}
	
	/**
     * Inner class for parameters
     */
    public static class VirialSWWEParam extends ParameterBase {
        public int nPoints = 5;
        public long numSteps = 1000;
        public static final int TREE = 0, CHAINS = 1, CHAIN_TAIL = 4, CHAIN_TREE = 5, CRINGS = 6, RING_TREE = 7, RINGS = 8, RING_CHAIN_TREES = 9;
        public int ref = RING_CHAIN_TREES;
        public double chainFrac = 0.3;
        public double treeFrac = 0.4;
        public double ringFrac = (1-chainFrac-treeFrac);
    }

}
