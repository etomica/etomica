/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.awt.Color;

import etomica.action.IAction;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.math.DoubleRange;
import etomica.data.histogram.HistogramSimple;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.ClusterWheatleySoft;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJ {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();
        if (args.length > 0) {
            params.writeRefPref = true;
        }
        ParseArgs.doParseArgs(params, args);
        if (args.length == 0) {
            params.nPoints = 4;
            params.temperature = 1;
            params.numSteps = 10000000L;
            params.doWheatley = true;
            params.doHist = false;
            params.doReeHoover = false;
        }
        
        runVirial(params);
    }

    public static void runVirial(VirialLJParam params) {
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        boolean doHist = params.doHist;
        boolean doWheatley = params.doWheatley;
        boolean doReeHoover = params.doReeHoover;

        final double HSBn = Standard.BHS(nPoints, sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSBn);
        System.out.println("Lennard Jones overlap sampling B"+nPoints+" at T="+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        P2LennardJones pTarget = new P2LennardJones(space);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        if (doWheatley) System.out.println("Wheatley");
        else System.out.println("Ree-Hoover: "+doReeHoover);
        VirialDiagrams diagrams = new VirialDiagrams(nPoints, false, false);
        diagrams.setDoReeHoover(true);
        diagrams.setDoShortcut(true);
        ClusterAbstract refCluster = doWheatley ? new ClusterWheatleyHS(nPoints, fRef) : diagrams.makeVirialCluster(fRef);
        refCluster.setTemperature(temperature);

        diagrams = new VirialDiagrams(nPoints, false, false);
        diagrams.setDoReeHoover(doReeHoover);
        diagrams.setDoShortcut(true);
        final ClusterAbstract targetCluster = doWheatley ? new ClusterWheatleySoft(nPoints, fTarget, 1e-12) : diagrams.makeVirialCluster(fTarget);
        targetCluster.setTemperature(temperature);

        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), temperature,refCluster,targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
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

        steps /= 1000;
        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        //sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/20);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/10);
        
        System.out.println("equilibration finished");
        
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
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

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

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

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(HSBn);
        System.out.println("time: "+(t2-t1)/1000.0);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 1;
        public long numSteps = 10000000;
        public double sigmaHSRef = 1.5;
        public boolean writeRefPref = false;
        public boolean doReeHoover = true;
        public double refFrac = -1;
        public boolean doHist = false;
        public boolean doWheatley = true;
    }
}
