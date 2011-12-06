package etomica.virial.simulations;


import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtomList;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IVectorMutable;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P2SphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.DoubleRange;
import etomica.util.HistogramNotSoSimple;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;

/**
 * Computes additive virial coefficients using the pair potential for He of Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.
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


public class VirialHeAJS {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; boolean QFH;
        double refFrac = params.refFrac;
        nPoints = params.nPoints;
        temperatureK = params.temperature;
        steps = params.numSteps;
        sigmaHSRef = params.sigmaHSRef;
        QFH = params.semiClassical;
            
        int numSubSteps = 1000;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);
        

        System.out.println("Target diagram: B"+nPoints+" for helium pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        if (QFH) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        System.out.println("B"+nPoints+"HS: "+HSB);
        //Next line not needed because energy in Kelvin
        //temperature = Kelvin.UNIT.toSim(temperature);

        
        Space space = Space3D.getInstance();

        MayerGeneralSpherical fTarget; MayerESpherical eTarget;
        Potential2Spherical p2 = new P2HePCKLJS(space);
        if (QFH) {
            P2EffectiveFeynmanHibbs p2QFH = new P2EffectiveFeynmanHibbs(Space3D.getInstance(),(Potential2SoftSpherical)p2);
			p2QFH.setTemperature(temperature);		
			p2QFH.setMass(4.002602);
			p2 = p2QFH;
        }
        if (params.truncation) {
            p2 = new P2SphericalTruncated(space, p2, 40);
        }
        fTarget = new MayerGeneralSpherical(p2);
        eTarget = new MayerESpherical(p2);
    	
    	ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
    	
   	    ClusterWeight sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster);

    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster);

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refSample.setTemperature(temperature);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), 
                temperature, new ClusterAbstract[]{refCluster,targetCluster},new ClusterWeight[]{refSample,sampleCluster1}, false);
        sim.integratorOS.setAgressiveAdjustStepFraction(true);

        IAtomList atoms = sim.box[1].getLeafList();
        double r = 4;
        for (int i=1; i<nPoints; i++) {
            IVectorMutable v = atoms.getAtom(i).getPosition();
            v.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
        
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
       
        sim.setAccumulatorBlockSize(steps);
        
        System.out.println();

        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesSpheresMono species = (SpeciesSpheresMono)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
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

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_2b_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/10); // 5000 IntegratorOverlap steps = 5e6 steps
        
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IIntegratorListener histListener = new IIntegratorListener() {
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
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
            
            public void integratorInitialized(IIntegratorEvent e) {}
        };
        FileWriter fw1 = null;
        if (true) {
            try {
                fw1 = new FileWriter("sample"+temperatureK+".dat");
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
            final FileWriter fw = fw1;
            IIntegratorListener sampleListener = new IIntegratorListener() {
                public void integratorStepStarted(IIntegratorEvent e) {}
                
                public void integratorStepFinished(IIntegratorEvent e) {
                    CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                    try {
                        fw.write(Math.sqrt(cPairs.getr2(0,1))+"\n");
                    }
                    catch (IOException ex) {
                        throw new RuntimeException(ex);
                    }
                }
                public void integratorInitialized(IIntegratorEvent e) {}
            };
            sim.integratorOS.getEventManager().addListener(sampleListener);
        }

        if (refFrac >= 0) {
            if (refFrac > 0.5) {
                sim.integrators[0].getEventManager().addListener(histListener);
            }
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }
        
        sim.integratorOS.getMoveManager().setEquilibrating(false);
//        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(steps);
        
        long t1 = System.currentTimeMillis();
        
        sim.getController().actionPerformed();
        if (true) {
            try {
                fw1.close();
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
        
        if (refFrac > 0.5) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            double[] hpi = piHist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i]) && h[i]!=0) {
                    System.out.println(xValues[i]+" "+h[i]+" "+hpi[i]);
                }
            }
        }
        
        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB);
        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

        public int nPoints = 2;
        public double temperature = 100;
        public long numSteps = 2000000000;
        public double sigmaHSRef = 3;
        public boolean semiClassical = false;
        public double refFrac = 0;
        public boolean truncation = false;
    }
}
