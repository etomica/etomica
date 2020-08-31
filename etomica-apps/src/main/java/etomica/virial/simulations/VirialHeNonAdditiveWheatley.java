/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;

import java.awt.*;
import java.util.Arrays;

/**
 * Adapted by Andrew from VirialHeNonAdditive
 * 
 * Computes only the nonadditive component of either the third, fourth, or fifth virial coefficient for the
 * ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009). 
 * 
 * Uses Wheatley recursion.
 */
public class VirialHeNonAdditiveWheatley {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize here
            params.calcApprox = true;
            params.nPoints = 6;
            params.numSteps = 1000000;
            params.semiClassical = true;
            params.subtractApprox = false;
            params.temperature = 273.15;
        }
        
    	final int nPoints = params.nPoints;
    	final double temperatureK = params.temperature;
        long steps = params.numSteps;
        final double sigmaHSRef = params.sigmaHSRef < 0 ? (3.0 + 120/(100+temperatureK)) : params.sigmaHSRef;
        final boolean semiClassical = params.semiClassical;
        final int nullRegionMethod = params.nullRegionMethod;
        double refFrac = params.refFrac;
        final boolean subtractApprox = params.subtractApprox;
        final boolean calcApprox = !subtractApprox && params.calcApprox;
        final double sigma = !calcApprox && !subtractApprox ? params.sigma : 0;
        final boolean doTotal = params.doTotal;
        
        double vhs = 4.0/3.0*Math.PI*Math.pow(sigmaHSRef, 3);
        final double HSBn = SpecialFunctions.factorial(nPoints)/2*Math.pow(vhs, nPoints-1);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSBn);
        System.out.println("Helium overlap sampling B"+nPoints+"NonAdd at T="+temperatureK+ " K");
        System.out.println("Using "+(semiClassical ? "semi" : "")+"classical pair potential");
        System.out.println("null region method = "+nullRegionMethod);
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (subtractApprox) {
            System.out.println("computing difference from approximate He");
        } else if (doTotal) {
            System.out.println("computing total");
        }

        final double temperature = Kelvin.UNIT.toSim(temperatureK);

        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        Space space = Space3D.getInstance();

        
        MayerFunction fRef = new MayerFunction() {

            public void setBox(Box box) {}
            public IPotential getPotential() {return null;}

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef*sigmaHSRef ? 1 : 0;
            }
        };

        MayerGeneralSpherical fTarget, fTargetSigma;
        MayerGeneralSpherical fTargetApprox;
        if (semiClassical) {
            P2HeSimplified p2cApprox = new P2HeSimplified(space);
            Potential2Spherical p2Approx = p2cApprox.makeQFH(temperature);
            
            P2HePCKLJS p2c = new P2HePCKLJS(space);
            Potential2Spherical p2 = p2c.makeQFH(temperature);
            P2HePCKLJS p2cs = new P2HePCKLJS(space, sigma);
            Potential2Spherical p2s = p2cs.makeQFH(temperature);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);
            fTargetSigma = new MayerGeneralSpherical(p2s);
        } else {
            P2HeSimplified p2Approx = new P2HeSimplified(space);
            
            P2HePCKLJS p2 = new P2HePCKLJS(space);
            P2HePCKLJS p2s = new P2HePCKLJS(space, sigma);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);
            fTargetSigma = new MayerGeneralSpherical(p2s);
        }

        IPotentialAtomicMultibody p3 = new P3CPSNonAdditiveHe(space);
        IPotentialAtomicMultibody p3sigma = new P3CPSNonAdditiveHe(space, sigma);
        P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);

        final MayerFunctionSphericalThreeBody f3Target = new MayerFunctionSphericalThreeBody(calcApprox ? p3Approx : p3);

        MayerFunctionNonAdditive[] fNA = new MayerFunctionNonAdditive[4];
        fNA[3] = f3Target;
        ClusterWheatleyMultibody fullTargetCluster = new ClusterWheatleyMultibody(nPoints, fTarget, fNA);
        fullTargetCluster.setTolerance(1e-14);
        fullTargetCluster.setDoTotal(doTotal);
        ClusterAbstract refCluster = new ClusterChainHS(nPoints, fRef);

        ClusterAbstract targetCluster = null;
        if (subtractApprox) {
            MayerFunctionSphericalThreeBody f3TargetApprox = new MayerFunctionSphericalThreeBody(p3Approx);
            MayerFunctionNonAdditive[] fNAapprox = new MayerFunctionNonAdditive[4];
            fNAapprox[3] = f3TargetApprox;
            final ClusterWheatleyMultibody[] targetSubtract = new ClusterWheatleyMultibody[1];
            targetSubtract[0] = new ClusterWheatleyMultibody(nPoints, fTargetApprox, fNAapprox);
            targetSubtract[0].setTolerance(1e-14);
            targetSubtract[0].setDoTotal(doTotal);
            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
        } else if (sigma != 0) {
            MayerFunctionSphericalThreeBody f3TargetSigma = new MayerFunctionSphericalThreeBody(p3sigma);
            MayerFunctionNonAdditive[] fNAsigma = new MayerFunctionNonAdditive[4];
            fNAsigma[3] = f3TargetSigma;
            ClusterWheatleyMultibody targetSigma = new ClusterWheatleyMultibody(nPoints, fTargetSigma, fNAsigma);
            targetSigma.setTolerance(1e-14);
            targetSigma.setDoTotal(doTotal);
            targetCluster = new ClusterDifference(targetSigma, new ClusterAbstract[]{fullTargetCluster});
        } else {
            targetCluster = fullTargetCluster;
        }

        targetCluster.setTemperature(temperature);
    	
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A"))), nPoints,
                temperature, refCluster,targetCluster);
        sim.init();
        int[] seeds = sim.getRandomSeeds();
        System.out.println("Random seeds: "+ Arrays.toString(seeds));


        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterAtomHSChain mcMoveHS = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHS);
        sim.accumulators[0].setBlockSize(1);

        ///////////////////////////////////////////////
        // Initialize non-overlapped configuration
        ///////////////////////////////////////////////
        
        IAtomList atoms = sim.box[1].getLeafList();
        double r = 3;
        for (int i=1; i<nPoints; i++) {
            Vector v = atoms.get(i).getPosition();
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
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesGeneral species = (SpeciesGeneral)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0, false);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
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
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }

        sim.setAccumulatorBlockSize((int)steps);
        sim.integratorOS.setNumSubSteps((int)steps);

        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        final HistogramSimple targHist = new HistogramSimple(200, new DoubleRange(-1, 6));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2 = cPairs.getr2(i, j);
                        double r = Math.sqrt(r2);
                        if (r > 1) {
                            r = Math.log(r);
                        }
                        else {
                            r -= 1;
                        }
                        targHist.addValue(r);
                    }
                }
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
                if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                if (Double.isInfinite(sim.dvo.getAverageAndError()[0])) {
                    sim.dvo.getAverageAndError();
                    throw new RuntimeException("oops");
                }
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSBn+", error: "+ratioAndError[1]*HSBn);
            }
        };
        if (!isCommandline) {
            sim.integratorOS.getEventManager().addListener(progressReport);
            if (params.doHist) {
                IntegratorListener histReport = new IntegratorListener() {
                    public void integratorInitialized(IntegratorEvent e) {}
                    public void integratorStepStarted(IntegratorEvent e) {}
                    public void integratorStepFinished(IntegratorEvent e) {
                        if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                        double[] xValues = targHist.xValues();
                        double[] h = targHist.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (!Double.isNaN(h[i]) && h[i] != 0) {
                                double r = xValues[i];
                                double y = h[i];
                                if (r < 0) r += 1;
                                else {
                                    r = Math.exp(r);
                                    y /= r;
                                }
                                System.out.println(r+" "+y);
                            }
                        }
                    }
                };
                sim.integratorOS.getEventManager().addListener(histReport);
            }
        }


        if (params.doHist) {
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }
        if (refFrac >= 0) {

            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }


        sim.integratorOS.getMoveManager().setEquilibrating(false);
sim.getController().runActivityBlocking(ai);
        
        long t2 = System.currentTimeMillis();
        
        if (params.doHist) {
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i]) && h[i] != 0) {
                    double rr = xValues[i];
                    double y = h[i];
                    if (rr < 0) rr += 1;
                    else {
                        rr = Math.exp(rr);
                        y /= rr;
                    }
                    System.out.println(rr+" "+y);
                }
            }
        }
        

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(HSBn);
        
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
        public double sigma = 0;
        public boolean doTotal = false;
    }
}
