/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.cluster.VirialDiagramsPT;

import javax.swing.*;
import java.awt.*;

public class VirialSwsPT {

    public static void main(String[] args) {
        IsomorphismFilter.DEBUG_MODE = false;
        VirialSiepmannSpheresParam params = new VirialSiepmannSpheresParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
        	// type here for running from Eclipse
        	params.nPoints = 4;
        	params.numSteps = 100000000;
        	params.lambda = 2;
        	params.order = 1;
        	params.doExp = true;
        }
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double sigmaHSRef = params.sigmaHSRef;
        int order = params.order;
        final double lambda = params.lambda;
        if (sigmaHSRef > lambda) {
            sigmaHSRef = lambda;
        }
        else if (sigmaHSRef <= 1.0) {
            sigmaHSRef = (lambda + 1.0) * 0.5;
        }
        boolean doExp = params.doExp;

        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fHS = new MayerHardSphere(sigmaHSRef);
        System.out.println("SW B"+nPoints+""+order);
        System.out.println("lambda: "+lambda);

        Potential2Spherical p2Ref = new P2HardSphere(space, 1.0,false);
    	Potential2Spherical p2Att = new P2SquareWell(space, 1.0, lambda, 1.0, false);
        
        MayerGeneralSpherical fTargetRef = new MayerGeneralSpherical(p2Ref);
        MayerSphericalPTAtt[] fTargetAtt = new MayerSphericalPTAtt[order];
        for (int i=0; i<fTargetAtt.length; i++) {
        	fTargetAtt[i] = new MayerSphericalPTAtt(p2Ref, p2Att, i+1);
        }

        VirialDiagramsPT alkaneDiagrams = new VirialDiagramsPT(nPoints, false, false);
        alkaneDiagrams.setDoExp(doExp);
        alkaneDiagrams.setOrderBeta(order);
        alkaneDiagrams.setDoShortcut(true);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = alkaneDiagrams.makeVirialCluster(fTargetRef, fTargetAtt);
        alkaneDiagrams = null;

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoShortcut(true);
        ClusterAbstract refCluster = rigidDiagrams.makeVirialCluster(fHS);

        final double refIntegral = HSB[nPoints];

        targetCluster.setTemperature(1.0);
        refCluster.setTemperature(1.0);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // overerr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        System.out.println(steps+" steps (1000 blocks of "+(steps/1000)+")");
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), ClusterWeightAbs.makeWeightCluster(targetCluster)};

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{species},
                new int[]{nPoints},1.0, new ClusterAbstract[]{refCluster, targetCluster}, sampleClusters, true);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

//        ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }
        steps /= 1000;

        ConfigurationClusterMove ccm = new ConfigurationClusterMove(space, sim.getRandom(), lambda);
        ccm.initializeCoordinates(sim.box[1]);
        if (sim.box[1].getSampleCluster().value(sim.box[1])==0) {
        	throw new RuntimeException("couldn't find an initial configuration");
        }

        sim.integratorOS.setNumSubSteps(1000);

        if(false) {
    double size = 10;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox1.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 1000, false);
    sim.equilibrate(null, 2000, false);
    sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B" + nPoints + " (L/mol)^" + (nPoints - 1));
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup, CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());

            IAction pushAnswer = new IAction() {
                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio * refIntegral;
                    averageBox.putData(data);
                    data.x = error * refIntegral;
                    errorBox.putData(data);
                }

                DataDouble data = new DataDouble();
            };
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints - 1}));
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
    return;
}
        
        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+""+order : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);

        
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);

        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes  "+sim.mcMoveTranslate[0].getStepSize()+" "+sim.mcMoveTranslate[1].getStepSize());

        if (false) {
            final double refIntegralF = refIntegral;
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    System.out.println("abs average: "+ratio*refIntegralF+", error: "+error*refIntegralF);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        long t0 = System.currentTimeMillis();
sim.getController().runActivityBlocking(ai);
        long t1 = System.currentTimeMillis();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(refIntegral);
        System.out.println("time: "+(t1-t0)*0.001);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialSiepmannSpheresParam extends ParameterBase {
        public int nPoints = 2;
        public long numSteps = 10000000;
        public double refFreq = -1;
        public int order = 1;
        public double lambda = 1.5;
        public double sigmaHSRef = 1.0;
        public boolean doExp = true;
    }
}
