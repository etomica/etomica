/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.MoleculeActionTranslateTo;
import etomica.api.ISpecies;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomList;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.IData;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.models.co2.PNGCPM;
import etomica.models.co2.PNGCPM.GCPMAgent;
import etomica.models.co2.PNGCPMX;
import etomica.models.water.SpeciesWater4PCOM;
import etomica.molecule.IMolecule;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.util.Arrays;

/**
 * Computes CO2-H2O mixture virial coefficients using ab-initio potentials
 * for both components.  Classical and semiclassical coefficients can be
 * computed.
 * 
 * 3-body dispersion for water (and water-CO2) can be used if Ewater is set.
 * CO2-H2O combining rules can be adjusted via k parameters.
 * 
 * @author Andrew Schultz
 */
public class VirialCO2H2OGCPMX {


    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.nTypes = new int[]{3,0};
            params.temperature = 280;
            params.numSteps = 10000000;
            params.sigmaHSRef = 5;
            params.nonAdditive = Nonadditive.FULL;
        }

    	final int[] nTypes = params.nTypes;
        final int nPoints = params.nTypes[0] + params.nTypes[1];
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;

        final double refFrac = params.refFrac;
        final Nonadditive nonAdditive = nPoints < 3 ? Nonadditive.NONE : params.nonAdditive;
        final double kijSigma = params.kijSigma;
        final double kijEpsilon = params.kijEpsilon;
        final double kijGamma = params.kijGamma;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);

        System.out.println("Overlap sampling for CO2/H2O GCPM mixture at " + temperatureK + " K");
        System.out.println("kEpsilon: "+kijEpsilon+"  kSigma: "+kijSigma+"  kGamma: "+kijGamma);
        if (nonAdditive != Nonadditive.NONE) {
            if (nonAdditive == Nonadditive.DISPERSION) {
                System.out.println("Including non-additive dispersion");
            }
            else {
                System.out.println("Including non-additive dispersion and induction");
            }
        }

        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB);
		
        final Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        
        SpeciesSpheresHetero speciesCO2 = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
        speciesCO2.setChildCount(new int[]{1,2});
        speciesCO2.setConformation(new IConformation() {
            
            public void initializePositions(IAtomList atomList) {
                atomList.getAtom(0).getPosition().E(0);
                atomList.getAtom(1).getPosition().setX(0,1.161);
                atomList.getAtom(2).getPosition().setX(0,-1.161);
            }
        });

        SpeciesWater4PCOM speciesWater = new SpeciesWater4PCOM(space);

        AtomTypeAgentManager paramsManager = new AtomTypeAgentManager(null);
        final PNGCPMX pTarget = new PNGCPMX(space, paramsManager, 6, kijSigma, kijEpsilon, kijGamma);

        MayerGeneral fTarget = new MayerGeneral(pTarget);

        ClusterAbstract targetCluster = new ClusterWheatleySoft(nPoints, fTarget, 1e-12);
        if (nTypes[1] == 2 && nPoints == 2) {
            // pure B2 for water.  we need flipping.
            // additive B3 for water should be fine and biconnectivity will help with mixture coefficients.
            ((ClusterWheatleySoft)targetCluster).setDoCaching(false);
            targetCluster = new ClusterCoupledFlipped(targetCluster, space, 20);
        }

        if (nonAdditive != Nonadditive.NONE) {
            IPotentialMolecular p3ATM = pTarget.makeAxilrodTeller();
            PotentialNonAdditive pi = null;
            IPotentialMolecular[] allPi = new IPotentialMolecular[nPoints-1];

            if (nonAdditive != Nonadditive.DISPERSION) {
                // we can only handle 3-body for now
                allPi[0] = pTarget.makeCachedPairPolarization();
                for (int i=1; i<allPi.length; i++) {
                    PNGCPM p = new PNGCPM(space, paramsManager, 6, 2+i);
                    p.setComponent(PNGCPM.Component.INDUCTION);
                    allPi[i] = p;
                }

                pi = new PotentialNonAdditive(allPi);

                MayerFunctionNonAdditiveFull[] fNA = new MayerFunctionNonAdditiveFull[nPoints+1];
                for (int i=3; i<=nPoints; i++) {
                    IPotentialMolecular p = null;
                    if (i==3 && nTypes[0] > 2) {
                        p = null;
                    }
                    else {
                        p = pi.makeNB(i);
                    }
                    fNA[i] = new MayerFunctionNonAdditiveFull(p);
                }
                targetCluster = new ClusterWheatleyMultibody(nPoints, fTarget, fNA);
            }
            else {
                MayerFunctionMolecularThreeBody f3 = new MayerFunctionMolecularThreeBody(p3ATM);
                targetCluster = new ClusterWheatleyMultibody(nPoints, fTarget, f3);
            }

            ((ClusterWheatleyMultibody)targetCluster).setRCut(100);
            if (nonAdditive == Nonadditive.FULL && nTypes[1] > 0) {
                // water induction requires flipping
                ((ClusterWheatleyMultibody)targetCluster).setDoCaching(false);
                targetCluster = new ClusterCoupledFlipped(targetCluster, space, 20);
            }
        }

        
        targetCluster.setTemperature(temperature);

        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
 		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{speciesCO2,speciesWater}, nTypes, temperature, refCluster, targetCluster);
//        sim.setRandom(new RandomMersenneTwister(new int[]{1941442288, -303985770, -1766960871, 2058398830}));
        sim.init();
        System.out.println("random seeds: "+Arrays.toString(sim.getRandomSeeds()));
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        paramsManager.setAgent(speciesWater.getHydrogenType(), new GCPMAgent(1.0,0,0.455,12.75,Electron.UNIT.toSim(0.6113),0,0,0));
        paramsManager.setAgent(speciesWater.getOxygenType(), new GCPMAgent(3.69,Kelvin.UNIT.toSim(110),0,12.75,0,0,0,0,0));
        paramsManager.setAgent(speciesWater.getMType(), new GCPMAgent(1.0,0,0.610,12.75,Electron.UNIT.toSim(-1.2226),0,0,0));
        paramsManager.setAgent(speciesWater.getCOMType(), new GCPMAgent(1.0,0,0.610,12.75,0,1.444,1.444,0));
        double qC = Electron.UNIT.toSim(0.6642);
        paramsManager.setAgent(speciesCO2.getAtomType(0), new GCPMAgent(3.193,Kelvin.UNIT.toSim(71.34),0.61/1.0483,15.5,qC,4.05,1.95,16.0/9.0*Kelvin.UNIT.toSim(2.52e4)) {
            protected final Vector r = space.makeVector();
            public Vector getParallelAxis(IMolecule mol) {
                IAtomList atoms = mol.getChildList();
                r.Ev1Mv2(atoms.getAtom(2).getPosition(),atoms.getAtom(1).getPosition());
                r.normalize();
                return r;
            }
        });
        double qO = -0.5*qC;
        paramsManager.setAgent(speciesCO2.getAtomType(1), new GCPMAgent(3.193*1.0483,Kelvin.UNIT.toSim(67.72),0.61,15.5,qO,0,0,0));
        
        if (nonAdditive != Nonadditive.NONE) {
            MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(space);
            Vector pos = space.makeVector();
            double r = 4;
            for (int i=1; i<nPoints; i++) {
                double theta = 2*i*Math.PI/nPoints;
                pos.setX(0, r*(1-Math.cos(theta)));
                pos.setX(1, r*Math.sin(theta));
                act.setDestination(pos);
                act.actionPerformed(sim.box[1].getMoleculeList().getMolecule(i));
            }
            sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
        }
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);

//            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
//            displayBox0.setColorScheme(colorScheme);
//            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
//            displayBox1.setColorScheme(colorScheme);
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
            refFileName = "refpref"+nPoints+"_"+(nonAdditive== Nonadditive.NONE ?"2":("3"+(nonAdditive==Nonadditive.DISPERSION?"a":"ai")))+"_"+tempString;
            refFileName += "C";
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        int nBins = 100;
        double dx = sigmaHSRef/nBins;
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IntegratorListener histListenerRef = new IntegratorListener() {
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
                double v = finalTargetCluster.value(sim.box[0]);
                hist.addValue(Math.sqrt(r2Max), v);
                piHist.addValue(Math.sqrt(r2Max), Math.abs(v));
            }
            
            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double v = finalTargetCluster.value(sim.box[1]);
                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                }
                else {
                    r -= 1;
                }
                targHist.addValue(r, v);
                targPiHist.addValue(r, Math.abs(v));
            }

            public void integratorInitialized(IntegratorEvent e) {}
        };

        if (params.doHist) {
            IntegratorListener histReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.println("**** reference ****");
                    double[] xValues = hist.xValues();
                    double[] h = hist.getHistogram();
                    double[] piH = piHist.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            System.out.println(xValues[i]+" "+h[i]+" "+piH[i]);
                        }
                    }
                    System.out.println("**** target ****");
                    xValues = targHist.xValues();
                    h = targHist.getHistogram();
                    piH = targPiHist.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            double r = xValues[i];
                            if (r < 0) r += 1;
                            else r = Math.exp(r);
                            System.out.println(r+" "+h[i]+" "+piH[i]);
                        }
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(histReport);

            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.getController().actionPerformed();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            
            System.out.println("final ref histogram");
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
//                    System.out.println(xValues[i]+" "+(-2*h[i]+1)+" "+Math.exp(-u/temperature));
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }


        System.out.println();
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println();
        
        sim.printResults(HSB);

        DataGroup allYourBase = (DataGroup)sim.accumulators[1].getData();
        IData averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        IData errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        IData covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        int n = 0;
        double correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("diagram "+nPoints+"bc average: %20.15e error: %9.4e ocor: %6.4f\n",
                averageData.getValue(0), errorData.getValue(0), correlationCoef));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
    }
    
    enum Nonadditive {
        NONE, DISPERSION, FULL
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public int[] nTypes = new int[]{2,1};
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public boolean doHist = false;
        public Nonadditive nonAdditive = Nonadditive.NONE;
        public double kijSigma = 0.99;
        public double kijEpsilon = 1.10;
        public double kijGamma = 0.96;
    }
}
