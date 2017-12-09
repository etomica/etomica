/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.IData;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import java.awt.*;
import java.util.Set;

/* 
 * Adapted by Kate from VirialGCPM
 * 
 * Computes only the nonadditive component of either the third, fourth, or fifth virial coefficient for the
 * ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009). 
 * 
 * 
 */


public class VirialHeNonAdditive {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize here
        }
        
    	final int nPoints = params.nPoints;
    	final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef < 0) {
            sigmaHSRef = 3.0 + 120/(100+temperatureK);
        }
        final boolean semiClassical = params.semiClassical;
        final int nullRegionMethod = params.nullRegionMethod;
        double refFrac = params.refFrac;
        final boolean subtractApprox = params.subtractApprox;
        final boolean calcApprox = !subtractApprox && params.calcApprox;
        final boolean minMulti = params.minMulti;
        
        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Helium overlap sampling B"+nPoints+"NonAdd at T="+temperatureK+ " K");
        System.out.println("Using "+(semiClassical ? "semi" : "")+"classical pair potential");
        System.out.println("null region method = "+nullRegionMethod);
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (subtractApprox) {
            System.out.println("computing difference from approximate He");
        }

        final double temperature = Kelvin.UNIT.toSim(temperatureK);

        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        Space space = Space3D.getInstance();

        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        
        MayerGeneralSpherical fTarget;
        MayerGeneralSpherical fTargetApprox;
        if (semiClassical) {
            P2HeSimplified p2cApprox = new P2HeSimplified(space);
            Potential2Spherical p2Approx = p2cApprox.makeQFH(temperature);
            
            P2HePCKLJS p2c = new P2HePCKLJS(space);
            Potential2Spherical p2 = p2c.makeQFH(temperature);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);

        } else {
            P2HeSimplified p2Approx = new P2HeSimplified(space);
            
            P2HePCKLJS p2 = new P2HePCKLJS(space);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);
        }

        IPotentialAtomicMultibody p3 = new P3CPSNonAdditiveHe(space);
        P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);

        final MayerFunctionSphericalThreeBody f3Target = new MayerFunctionSphericalThreeBody(calcApprox ? p3Approx : p3);
        

        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, false);
        flexDiagrams.setDoMinimalMulti(true);
        flexDiagrams.setDoMultiFromPair(!minMulti);
        flexDiagrams.setDoReeHoover(true);
        flexDiagrams.setDoShortcut(true);
        ClusterSum fullTargetCluster = flexDiagrams.makeVirialCluster(fTarget, f3Target, false);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        rigidDiagrams.setAllPermutations(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);


        ClusterAbstract targetCluster = null;
        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];
        if (subtractApprox) {
            final ClusterSum[] targetSubtract = new ClusterSum[1];
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            MayerFunctionSphericalThreeBody f3TargetApprox = new MayerFunctionSphericalThreeBody(p3Approx);
            targetSubtract[0] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetApprox}, new MayerFunctionNonAdditive[]{f3TargetApprox});
            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
            ClusterSum[] targetDiagramsPlus = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)fullTargetCluster, fTarget, f3Target);
            targetDiagrams = new ClusterDifference[targetDiagramsPlus.length];
            ClusterSum[][] targetDiagramsMinus = new ClusterSum[targetDiagramsPlus.length][0];
            for (int j=0; j<targetDiagramsMinus.length; j++) {
                targetDiagramsMinus[j] = new ClusterSum[targetSubtract.length];
            }
            ClusterSum[] foo = null;
            foo = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)targetSubtract[0], fTargetApprox, f3TargetApprox);
            for (int j=0; j<foo.length; j++) {
                targetDiagramsMinus[j][0] = foo[j];
            }
            for (int j=0; j<targetDiagramsPlus.length; j++) {
                targetDiagrams[j] = new ClusterDifference(targetDiagramsPlus[j], targetDiagramsMinus[j]);
            }
        }
        else {
            targetCluster = fullTargetCluster;
            targetDiagrams = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)fullTargetCluster, fTarget, f3Target);
        }

        targetCluster.setTemperature(temperature);
    	
        refCluster.setTemperature(temperature);

        int[] targetDiagramNumbers = new int[targetDiagrams.length];
        int[] mfTargetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = flexDiagrams.getMSMCGraphs(true, true);
        int iGraph = 0;
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ed = new DeleteEdgeParameters(flexDiagrams.mmBond);
        for (Graph g : singleGraphs) {
            if (g.edgeCount() == nPoints*(nPoints-1)/2) {
                if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.mmBond)) {
                    System.out.print(" ("+g.coefficient()+") "+g.nodeCount()+"M");
                    targetDiagramNumbers[iGraph] = -g.nodeCount();
                }
            }
            else {
                String gnStr = g.getStore().toNumberString();
                targetDiagramNumbers[iGraph] = Integer.parseInt(gnStr);
                Graph gOnlyF = edgeDeleter.apply(g, ed);
                gnStr += "m"+gOnlyF.getStore().toNumberString();
                mfTargetDiagramNumbers[iGraph] = Integer.parseInt(gOnlyF.getStore().toNumberString());
                System.out.print(" ("+g.coefficient()+") "+gnStr);
            }
            System.out.println();
            iGraph++;
        }
        System.out.println();
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }
        ClusterWeight targetSampleCluster = ClusterWeightAbs.makeWeightCluster(targetCluster);
        ClusterWeight refSampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);



        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{new SpeciesSpheresMono(space, new ElementSimple("A"))}, new int[]{nPoints},
                temperature, new ClusterAbstract[]{refCluster,targetCluster}, targetDiagrams, new ClusterWeight[]{refSampleCluster,targetSampleCluster}, false);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        ///////////////////////////////////////////////
        // Initialize non-overlapped configuration
        ///////////////////////////////////////////////
        
        IAtomList atoms = sim.box[1].getLeafList();
        double r = 4;
        for (int i=1; i<nPoints; i++) {
            Vector v = atoms.getAtom(i).getPosition();
            v.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
        
        /*
        IAtomList atoms0 = sim.box[0].getLeafList();
        for (int i=1;i<atoms0.getAtomCount();i++) {
        	atoms0.getAtom(i).getPosition().setX(0, i*sigmaHSRef*1.3);
        } */ 
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesSpheresMono species = (SpeciesSpheresMono)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }

        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = null;
        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints+"_3b_"+tempString;
            refFileName += semiClassical ? "_sc" : "_c";
            if (calcApprox) {
                refFileName += "a";
            }
            else if (subtractApprox) {
                refFileName += "sa";
            }
        }
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);

        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize((int)steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());
        
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IntegratorListener histListener = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                CoordinatePairSet cPairs = sim.box[0].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }
                double v = sim.box[1].getSampleCluster().value(sim.box[0]);
                hist.addValue(Math.sqrt(r2Max), v);
                piHist.addValue(Math.sqrt(r2Max), Math.abs(v));
            }
            
            public void integratorInitialized(IntegratorEvent e) {}
        };
        IntegratorListener progressReport = new IntegratorListener() {
            public void integratorInitialized(IntegratorEvent e) {}
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorStepFinished(IntegratorEvent e) {
//                if (Double.isInfinite(sim.dsvo.getOverlapAverageAndError()[0])) {
//                    sim.dsvo.getOverlapAverageAndError();
//                    throw new RuntimeException("oops");
//                }
                if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                if (Double.isInfinite(sim.dvo.getAverageAndError()[0])) {
                    sim.dvo.getAverageAndError();
                    throw new RuntimeException("oops");
                }
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
            }
        };
        if (!isCommandline) {
            sim.integratorOS.getEventManager().addListener(progressReport);
            if (params.doHist) {
                IntegratorListener histReport = new IntegratorListener() {
                    public void integratorInitialized(IntegratorEvent e) {}
                    public void integratorStepStarted(IntegratorEvent e) {}
                    public void integratorStepFinished(IntegratorEvent e) {
                        if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                        double[] xValues = hist.xValues();
                        double[] h = hist.getHistogram();
                        double[] piH = piHist.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (!Double.isNaN(h[i])) {
                                System.out.println(xValues[i]+" "+h[i]+" "+piH[i]);
                            }
                        }
                    }
                };
                sim.integratorOS.getEventManager().addListener(histReport);
            }
        }
        

        if (refFrac >= 0) {
            if (params.doHist) {
                sim.integrators[0].getEventManager().addListener(histListener);
            }
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }


        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(1000);
        sim.getController().actionPerformed();
        
        long t2 = System.currentTimeMillis();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i]) && h[i]!=0) {
                    System.out.println(xValues[i]+" "+h[i]);
                }
            }
        }
        

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(HSB[nPoints]);

        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(AccumulatorAverage.AVERAGE.index);
        IData dataErr = allData.getData(AccumulatorAverage.ERROR.index);
        IData dataCov = allData.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        int nTotal = (targetDiagrams.length+2);
        double oVar = dataCov.getValue(nTotal*nTotal-1);
        for (int i=0; i<targetDiagrams.length; i++) {
            if (targetDiagramNumbers[i]<0) {
                System.out.print("diagram "+(-targetDiagramNumbers[i])+"M ");
            }
            else {
                if (mfTargetDiagramNumbers[i]<0) {
                    System.out.print("diagram "+targetDiagramNumbers[i]+"rm"+(-mfTargetDiagramNumbers[i])+" ");
                }
                else {
                    System.out.print("diagram "+targetDiagramNumbers[i]+"m"+mfTargetDiagramNumbers[i]+" ");
                }
            }
            // average is vi/|v| average, error is the uncertainty on that average
            // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)
            double ivar = dataCov.getValue((i+1)*nTotal+(i+1));
            double ocor = ivar*oVar == 0 ? 0 : dataCov.getValue(nTotal*(i+1)+nTotal-1)/Math.sqrt(ivar*oVar);
            System.out.print(String.format("average: %20.15e  error: %10.15e  ocor: %7.5f", dataAvg.getValue(i+1), dataErr.getValue(i+1), ocor));
            if (targetDiagrams.length > 1) {
                System.out.print("  dcor:");
                for (int j=0; j<targetDiagrams.length; j++) {
                    if (i==j) continue;
                    double jvar = dataCov.getValue((j+1)*nTotal+(j+1));
                    double dcor = ivar*jvar == 0 ? 0 : dataCov.getValue((i+1)*nTotal+(j+1))/Math.sqrt(ivar*jvar);
                    System.out.print(String.format(" %8.6f", dcor));
                }
            }
            System.out.println();
        }
        
        System.out.println();
        System.out.println("time: "+(t2-t1)/1000.0);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these!
        public int nPoints = 4;
        public double temperature = 300;   // Kelvin
        public long numSteps = 1000000;
        public double sigmaHSRef = -1; // negative means use equation
        public int nullRegionMethod = 2;
        public double refFrac = -1;
        public boolean doHist = false;
        public boolean semiClassical = false;
        public boolean calcApprox = false;
        public boolean subtractApprox = false;
        public boolean minMulti = false;
    }
}
