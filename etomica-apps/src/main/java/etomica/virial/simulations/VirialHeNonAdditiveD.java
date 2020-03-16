/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
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
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.util.Arrays;

/**
 * Adapted by Andrew from VirialHeNonAdditive
 * 
 * Computes only the nonadditive component of either the third, fourth, or fifth virial coefficient for the
 * ab initio non-additive trimer potential1 for He developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009).
 * 
 * Uses Wheatley recursion.
 */
public class VirialHeNonAdditiveD {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize here
            params.nPoints = 4;
            params.nDer = 3;
            params.temperature = 15;
            params.numSteps = 1000000;
            params.potential = PotentialChoice.OLD;
            params.calcDiff = PotentialChoice.NONE;
            params.doTotal=false;
            params.dop3Aprox = false;
            params.dop3AproxDiff = false;
            params.BDtol = 1e-12; //1e-14 used previously
            params.seed = null;

            params.refFrac = -1;
            params.sigmaHSRef = 5;
            params.dorefpref = false;

            params.doHist = false;
            params.doChainRef = true;
        }

        final int nPoints = params.nPoints;
        final int nDer = params.nDer;
        double temperatureK = params.temperature;
        long steps = params.numSteps;
        final PotentialChoice potential = params.potential;
        final PotentialChoice calcDiff = params.calcDiff;
        final boolean doTotal = params.doTotal;
        final boolean dop3Aprox = params.dop3Aprox;
        final boolean dop3AproxDiff = params.dop3AproxDiff;
        final double BDtol = params.BDtol;
        int[] seed = params.seed;

        double refFrac = params.refFrac;
        final double sigmaHSRef = params.sigmaHSRef;
        boolean dorefpref = params.dorefpref;

        boolean doHist = params.doHist;
        boolean doChainRef = params.doChainRef;

        final double temperature = Kelvin.UNIT.toSim(temperatureK);

        long blockSize = 1000;
        int EqSubSteps = 1000;


        System.out.println("Overlap sampling for He at " + temperatureK + " K " + "for Nonadditive B"+nPoints+" and "+nDer+" derivatives");

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);

        String coeftype = doTotal ? "total" : "non-additive";

        if ( potential == PotentialChoice.NONE ) {
            throw new RuntimeException("POTENTIAL NOT FOUND !!!!");
        }
        else if ( calcDiff == PotentialChoice.NEW ){
            throw new RuntimeException("CALCDIFF CANNOT BE NEW !!!!");
        }
        else if ( potential == calcDiff && dop3Aprox == dop3AproxDiff) {
            throw new RuntimeException("POTENTIAL AND CALCDIFF ARE THE SAME !!!!");
        }
        else if ( calcDiff == PotentialChoice.NONE ){
            System.out.println("Computing full "+ coeftype +" coefficient using " + potential + " Potential");
            if(dop3Aprox){
                System.out.println("Using approximate 3-Body potential in main cluster");
            }
            else{
                System.out.println("Using rigorous 3-Body potential in main cluster");
            }
            if ( dop3AproxDiff ) {
                throw new RuntimeException(" NO DIFFERENCE TO APPLY P3 APPROX !!!! SET dop3AproxDiff TO FALSE OR SET CALCDIFF TO SIMPLE/OLD/NEW");
            }
        }
        else if ( calcDiff != PotentialChoice.NONE ) {
            if ( potential == PotentialChoice.SIMPLE && calcDiff == PotentialChoice.OLD ){
                throw  new RuntimeException("SWITCH POTENTIAL AND CALCDIFF !!!!");
            }
            else {
                System.out.println("Computing "+ coeftype + " coefficient difference between " + potential + " and " + calcDiff + " Potential");
            }
            if(dop3Aprox){
                System.out.println("Using approximate 3-Body potential in main cluster");
            }
            else{
                System.out.println("Using rigorous 3-Body potential in main cluster");
            }
            if(dop3AproxDiff){
                System.out.println("Using approximate 3-Body potential in difference cluster");
            }
            else{
                System.out.println("Using rigorous 3-Body potential in difference cluster");
            }
        }

        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        if (doChainRef) System.out.println("HS Chain reference");
        System.out.println("  B"+nPoints+"HS: "+HSBn);

        Space space = Space3D.getInstance();

        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public IPotential getPotential() {
                return null;
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        Potential2SoftSpherical pTarget ;

        if ( potential == PotentialChoice.SIMPLE ) {
            pTarget = new P2HeSimplified(space);
        }
        else if ( potential==PotentialChoice.OLD ) {
            pTarget = new P2HePCKLJS(space);
        }
        else {
            pTarget = new P2HePCJS(space);
        }

        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);

        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        IPotentialAtomicMultibody p3 = new P3CPSNonAdditiveHe(space);
        P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);

        final MayerFunctionSphericalThreeBody f3Target = new MayerFunctionSphericalThreeBody(dop3Aprox ? p3Approx : p3);
        
        MayerFunctionNonAdditive[] fNA = new MayerFunctionNonAdditive[4];
        fNA[3] = f3Target;
        final ClusterWheatleyMultibodyDerivatives fullTargetCluster = new ClusterWheatleyMultibodyDerivatives(nPoints, fTarget, fNA, BDtol, nDer,doTotal);

        ClusterAbstract targetCluster = null;
        ClusterWheatleyMultibodyDerivatives clusterDiff = null;
        if (calcDiff != PotentialChoice.NONE) {
            Potential2SoftSpherical pTargetDiff;
            if ( calcDiff == PotentialChoice.OLD ) {
                pTargetDiff = new P2HePCKLJS(space);
            }
            else {
                pTargetDiff = new P2HeSimplified(space);
            }
            MayerGeneralSpherical fTargetDiff = new MayerGeneralSpherical(pTargetDiff);
            MayerFunctionSphericalThreeBody f3TargetDiff = new MayerFunctionSphericalThreeBody(dop3AproxDiff ? p3Approx : p3);
            MayerFunctionNonAdditive[] fNADiff = new MayerFunctionNonAdditive[4];
            fNADiff[3] = f3TargetDiff;
            ClusterAbstract[] targetSubtract = new ClusterAbstract[1];
            clusterDiff = new ClusterWheatleyMultibodyDerivatives(nPoints, fTargetDiff, fNADiff, BDtol, nDer,doTotal);;
            targetSubtract[0] = clusterDiff;
            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
        }
        else {
            targetCluster = fullTargetCluster;
        }

        targetCluster.setTemperature(temperature);

        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), nPoints, temperature,refCluster,targetCluster);
        if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));

        ClusterAbstract[] targetDiagrams = null;
        if(calcDiff != PotentialChoice.NONE){
            targetDiagrams = new ClusterDifference[nDer];
            for(int m=1;m<=nDer;m++){
                ClusterAbstract targetClusterm = new ClusterWheatleyMultibodyDerivatives.ClusterRetrievePrimes(fullTargetCluster,m);
                ClusterAbstract[] targetClusterDiffm = new ClusterAbstract[1];
                targetClusterDiffm[0] = new ClusterWheatleyMultibodyDerivatives.ClusterRetrievePrimes(clusterDiff,m);
                targetDiagrams[m-1]= new ClusterDifference(targetClusterm,targetClusterDiffm);
            }
        }
        else{
            targetDiagrams = new ClusterWheatleyMultibodyDerivatives.ClusterRetrievePrimes[nDer];
            for(int m=1;m<=nDer;m++){
                targetDiagrams[m-1]= new ClusterWheatleyMultibodyDerivatives.ClusterRetrievePrimes(fullTargetCluster,m);
            }
        }

        sim.setExtraTargetClusters(targetDiagrams);
        sim.init();

        System.out.println("random seeds: "+ Arrays.toString(seed==null?sim.getRandomSeeds():seed));

        if (doChainRef) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.accumulators[0].setBlockSize(1);
        }


        ///////////////////////////////////////////////
        // Initialize non-overlapped configuration
        ///////////////////////////////////////////////

        IAtomList atoms = sim.box[1].getLeafList();
        double rt = 4;
        for (int i=1; i<nPoints; i++) {
            Vector v = atoms.get(i).getPosition();
            v.setX(0, rt*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, rt*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();

        ///////////////////////////////////////////////
        ////////////////////DONE///////////////////////
        ///////////////////////////////////////////////

        sim.integratorOS.setNumSubSteps(EqSubSteps);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            return;
        }

        long t1 = System.currentTimeMillis();

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        System.out.println();
        String refFileName = null;

        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            String coefStr = doTotal ? "_tot_" : "_nad_";
            String potStr = potential==PotentialChoice.SIMPLE?"simple_":(potential==PotentialChoice.OLD?"old_":(potential==PotentialChoice.NEW?"new_": " ???? ")) ;
            String cdStr = calcDiff==PotentialChoice.SIMPLE?"simple_diff_":(calcDiff==PotentialChoice.OLD?"old_diff_":(calcDiff==PotentialChoice.NONE?"full_": " ???? ")) ;
            String p3AStr = dop3Aprox?"p3A_":"";
            String p3ADStr = dop3AproxDiff?"p3AD_":"";
            refFileName = "refpref_"+nPoints+coefStr+potStr+cdStr+p3AStr+p3ADStr+tempString;
            refFileName += "K";
        }


        final HistogramSimple targHist = new HistogramSimple(200, new DoubleRange(-1, 4));
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

        if (doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.initRefPref(refFileName, (steps / EqSubSteps) / 20);
        sim.equilibrate(refFileName, (steps / EqSubSteps) / 10);

        System.out.println("equilibration finished");

        if(dorefpref){
            long t2 = System.currentTimeMillis();
            System.out.println("time: "+(t2-t1)/1000.0);
            return;
        }

        IntegratorListener progressReport = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                if (sim.integratorOS.getStepCount() % 100 != 0) return;
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSBn+", error: "+ratioAndError[1]*HSBn);
            }

            public void integratorInitialized(IntegratorEvent e) {}

        };
        if (false) {
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);

        if (doChainRef) sim.accumulators[0].setBlockSize(1);
        sim.ai.setMaxSteps(steps / blockSize);
        for (int i=0; i<2; i++) {
            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }

        sim.getController().actionPerformed();

        if (doHist) {
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
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

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());

        String[] extraNames = new String[nDer];
        for (int i = 1; i <= nDer; i++) {
            extraNames[i - 1] = "derivative " + i;
        }
        sim.printResults(HSBn, extraNames);

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
    }

    enum PotentialChoice {
        NONE, SIMPLE,  OLD, NEW
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these!

        public int nPoints = 2;
        public int nDer = 2;
        public double temperature = 250;          // Kelvin
        public long numSteps = 10000000;
        public PotentialChoice potential = PotentialChoice.NEW;
        public PotentialChoice calcDiff = PotentialChoice.NONE;
        public boolean doTotal = false;
        public boolean dop3Aprox = false;
        public boolean dop3AproxDiff = false;
        public double BDtol = 1e-12;
        public int[] seed = null;

        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public boolean dorefpref = false;

        public boolean doHist = false;
        public boolean doChainRef = true;
    }
}
