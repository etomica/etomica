/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.histogram.HistogramExpanding;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.starpolymer.ConformationStarPolymerGraft;
import etomica.units.Pixel;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.integrator.IntegratorMC;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Mayer sampling for star/grafted particles (baseline: star polymer).
 * Minimal changes:
 *  - distinct core vs monomer atom types
 *  - tunable core size (sigma0) fed to ConformationStarPolymerGraft
 *  - separate LJ params for core-core, mono-mono, core-mono
 *  - HS ref diameter optionally includes core size
 */
public class Virialjack {

    public static ClusterSum makeB2Cluster(MayerFunction f) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

            params.nPoints = 2;
            params.armLength =128;
            params.coreSigma=4;
            params.numArms = 20;
            params.temperature = 5;
            params.numSteps = 100000L;
            params.rc=2.5;
        }

        final int nPoints = params.nPoints;
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        if (armLength< 2){
            throw new RuntimeException("arm length must be at least 2");
        }

        Space space = Space3D.getInstance();

        // --- Conformation with tunable core size (sigma0) ---
        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(Space3D.getInstance(), numArms, armLength);

        conf.setSigma0(params.coreSigma);


        // --- Two atom types: core vs monomer (keeps indices identical) ---
        AtomType typeCore = AtomType.simple("Core");
        AtomType typeMono = AtomType.simple("Mono");

        // species: 1 core + (numArms*armLength) monomers
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeCore, 1)
                .addCount(typeMono, numArms * armLength)
                .withConformation(conf)
                .build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        // Pair potentials (molecule-level wrapper)
        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);

        System.out.println(numArms+" arms of length "+armLength+"  B"+nPoints+" at T = "+temperature);

        // --- Define LJ potentials per pair (defaults preserve old behavior: ε=1, σ=1) ---
        // --- Define LJ using combining rules (Lorentz–Berthelot) ---
        double sigC = params.coreSigma;   // core σ
        double epsC = params.epsCore;     // core ε
        double sigM = 1;   // mono σ
        double epsM = 1;     // mono ε

        double sigCM = 0.5*(sigC + sigM);           // Lorentz
        double epsCM = Math.sqrt(epsC * epsM);      // Berthelot

        TruncationFactory tf = params.rc==0?new TruncationFactoryNone():new TruncationFactorySimple(params.rc);
        IPotential2 p2CoreCore  = P2LennardJones.makeTruncated(sigC,epsC,tf);
        IPotential2 p2MonoMono = P2LennardJones.makeTruncated(sigM,epsM,tf);
        IPotential2 p2CoreMono = P2LennardJones.makeTruncated(sigCM,epsCM,tf);


          // Mayer functions
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        // HS reference diameter: original heuristic plus optional core contribution
        // original: 1.5 + 0.15*(2*armLength)
        double sigmaHSRef = params.useCoreInHSRef
                ? (params.coreSigma + 0.15*(2*armLength) + params.hsExtra)
                : (1.5 + 0.15*(2*armLength) + params.hsExtra);

        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {}
            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        boolean flex = nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, flex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        // --- Bonding info (unchanged logic; just ensure intramolecular LJ is computed) ---
        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) { return false; }
        };
        // dummy bonding potential to force intramolecular nonbonded calc
        IPotential2 pBonding = new IPotential2() {
            @Override public double getRange() { return 2; }
            @Override public void u012add(double r2, double[] u012) { }
        };
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<numArms; i++) {
            pairs.add(new int[]{0,1+armLength*i});
            for (int j = 0; j < armLength-1; j++) {
                pairs.add(new int[]{1+armLength*i+j, 1+armLength*i+j+1});
            }
        }
        bondingInfo.setBondingPotentialPair(species, pBonding, pairs);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm,
                new int[]{flex ? (nPoints+1) : nPoints}, temperature, refCluster, targetCluster);
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.setRandom(new RandomMersenneTwister(3));
        sim.init();
        pTarget.integrator1=sim.integrators[0];
        pTarget.integrator2=sim.integrators[1];
        System.out.println("random seeds: "+ Arrays.toString(sim.getRandomSeeds()));

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
        sim.accumulators[0].setBlockSize(1);

        int[] constraintMap = new int[nPoints+1];
        for (int i=0; i<nPoints; i++) { constraintMap[i] = i; }
        if (flex) {
            constraintMap[nPoints] = 0;
            mcMoveHSC.setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
        }

        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        sim.integratorOS.setNumSubSteps(1000);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        // --- Assign pair potentials to atom-type pairs (key change) ---
        // Default old behavior recovered if eps/sig all = 1.
        pTarget.setAtomPotential(typeCore, typeCore, p2CoreCore);
        pTarget.setAtomPotential(typeMono, typeMono, p2MonoMono);
        pTarget.setAtomPotential(typeCore, typeMono, p2CoreMono);
        pTarget.setAtomPotential(typeMono, typeCore, p2CoreMono);

        sim.integratorOS.setNumSubSteps(1000);

        MCMoveClusterAngle[] angleMoves;
        MCMoveClusterStretch[] stretchMoves = null;

        IntArrayList[] bonding = new IntArrayList[1+numArms*armLength];
        bonding[0] = new IntArrayList(numArms);
        for (int i=0; i<numArms; i++) {
            bonding[0].add(1+i*armLength);
            if (armLength > 1) {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0, 1 + i * armLength + 1});
            } else {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0});
            }
            for (int j = 1; j < armLength-1; j++) {
                bonding[1+armLength*i+j] = new IntArrayList(new int []{1+armLength*i+j-1, 1+armLength*i+j+1});
            }
            if (armLength > 1) {
                bonding[1 + (i+1) * armLength - 1] = new IntArrayList(new int[]{1 + (i+1) * armLength - 2});
            }
        }

        PotentialCompute pc0 = sim.integrators[0].getPotentialCompute();
        PotentialCompute pc1 = sim.integrators[1].getPotentialCompute();

        /*angleMoves = new MCMoveClusterAngle[2];
        angleMoves[0] = new MCMoveClusterAngle(pc0, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength));
        angleMoves[0].setBox(sim.box[0]);
        angleMoves[0].setConstraintMap(constraintMap);
        sim.integrators[0].getMoveManager().addMCMove(angleMoves[0]);
        angleMoves[1] = new MCMoveClusterAngle(pc1, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength));
        angleMoves[1].setBox(sim.box[1]);
        angleMoves[1].setConstraintMap(constraintMap);
        sim.integrators[1].getMoveManager().addMCMove(angleMoves[1]);*/

        angleMoves = new MCMoveClusterAngle[4];
        angleMoves[0] = new MCMoveClusterAngle(pc0, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,1,armLength/2));
        angleMoves[0].setBox(sim.box[0]);
        angleMoves[0].setConstraintMap(constraintMap);
        angleMoves[0].setFixedCOM(false);
sim.integrators[0].getMoveManager().addMCMove(angleMoves[0]);
        angleMoves[1] = new MCMoveClusterAngle(pc1, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,1,armLength/2));
        angleMoves[1].setBox(sim.box[1]);
        angleMoves[1].setConstraintMap(constraintMap);
        angleMoves[1].setFixedCOM(false);
sim.integrators[1].getMoveManager().addMCMove(angleMoves[1]);

        angleMoves[2] = new MCMoveClusterAngle(pc0, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,armLength/2+1,armLength));
        angleMoves[2].setBox(sim.box[0]);
        angleMoves[2].setConstraintMap(constraintMap);
        angleMoves[2].setFixedCOM(false);
sim.integrators[0].getMoveManager().addMCMove(angleMoves[2]);
        angleMoves[3] = new MCMoveClusterAngle(pc1, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,armLength/2+1,armLength));
        angleMoves[3].setBox(sim.box[1]);
        angleMoves[3].setConstraintMap(constraintMap);
        angleMoves[3].setFixedCOM(false);
sim.integrators[1].getMoveManager().addMCMove(angleMoves[3]);

        MCMoveClusterAngleMulti angleMovesMulti1 = new MCMoveClusterAngleMulti(pc1, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,1,armLength/2),5);
        angleMovesMulti1.setBox(sim.box[1]);
        angleMovesMulti1.setConstraintMap(constraintMap);
        angleMovesMulti1.setFixedCOM(false);
        sim.integrators[1].getMoveManager().addMCMove(angleMovesMulti1);
        //new
        MCMoveClusterAngleMulti angleMovesMulti0 = new MCMoveClusterAngleMulti(pc0, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,1,armLength/2),5);
        angleMovesMulti0.setBox(sim.box[0]);
        angleMovesMulti0.setConstraintMap(constraintMap);
        angleMovesMulti0.setFixedCOM(false);
        sim.integrators[0].getMoveManager().addMCMove(angleMovesMulti0);

        MCMoveClusterAngleMulti angleMovesMulti0outer = new MCMoveClusterAngleMulti(pc0, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,armLength/2+1,armLength),10);
        angleMovesMulti0outer.setBox(sim.box[0]);
        angleMovesMulti0outer.setConstraintMap(constraintMap);
        angleMovesMulti0outer.setFixedCOM(false);
        sim.integrators[0].getMoveManager().addMCMove(angleMovesMulti0outer);//tillthis

        MCMoveClusterAngleMulti angleMovesMulti1outer = new MCMoveClusterAngleMulti(pc1, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(),numArms,armLength,armLength/2 + 1,armLength),10);
        angleMovesMulti1outer.setBox(sim.box[1]);
        angleMovesMulti1outer.setConstraintMap(constraintMap);
        angleMovesMulti1outer.setFixedCOM(false);
        sim.integrators[1].getMoveManager().addMCMove(angleMovesMulti1outer);


        MCMoveClusterShuffle shuffleMove0 = null, shuffleMove1 = null;
        if (false && armLength > 5 ) {
            shuffleMove0 = new MCMoveClusterShuffle(pc0, space, sim.getRandom());
            shuffleMove0.setBox(sim.box[0]);
            shuffleMove0.setBonding(bonding);
            shuffleMove0.setStepSizeMax(armLength - 2);
            shuffleMove0.setConstraintMap(constraintMap);
            //sim.integrators[0].getMoveManager().addMCMove(shuffleMove0);
            ((MCMoveStepTracker) shuffleMove0.getTracker()).setAcceptanceTarget(0.3);
            shuffleMove1 = new MCMoveClusterShuffle(pc1, space, sim.getRandom());
            shuffleMove1.setBonding(bonding);
            shuffleMove1.setBox(sim.box[1]);
            shuffleMove1.setStepSizeMax(armLength - 2);
            shuffleMove1.setConstraintMap(constraintMap);
            //sim.integrators[1].getMoveManager().addMCMove(shuffleMove1);
            ((MCMoveStepTracker) shuffleMove1.getTracker()).setAcceptanceTarget(0.3);
        }

        // (interactive visualization block unchanged and disabled)

        /*long t1 = System.nanoTime();

        targetCluster.setTemperature(temperature * 100);
        long relaxSteps = steps/40;
        Activity activityInitRefPref = new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                angleMoves[0].skipW=true;
                angleMoves[1].skipW=true;
                for (int i = 0; i < relaxSteps*1000; i++) {

                    angleMoves[0].doTrial();
                    double chi=angleMoves[0].getChi(temperature*100);
                    if(chi<sim.getRandom().nextDouble()){
                        angleMoves[0].rejectNotify();
                    }
                    else {
                        angleMoves[0].acceptNotify();
                    }
                    angleMoves[1].doTrial();
                    chi=angleMoves[1].getChi(temperature*100);
                    if(chi<sim.getRandom().nextDouble()){
                        angleMoves[1].rejectNotify();
                    }
                    else {
                        angleMoves[1].acceptNotify();
                    }

                }
                angleMoves[0].skipW=false;
                angleMoves[1].skipW=false;
            }
        };
        sim.getController().runActivityBlocking(activityInitRefPref);
        targetCluster.setTemperature(temperature);*/

        long t1 = System.nanoTime();

        targetCluster.setTemperature(temperature * 100);
        long relaxSteps = steps/40;
        IntegratorMC relaxIntegrator= new IntegratorMC(sim.integrators[1].getPotentialCompute(), sim.getRandom(), temperature*100,sim.box[1]);
        angleMoves[1].skipW=true;
        angleMoves[3].skipW=true;
        angleMovesMulti1.skipW=true;
        angleMovesMulti1outer.skipW=true;
        relaxIntegrator.getMoveManager().addMCMove(angleMoves[1],0.1);
        relaxIntegrator.getMoveManager().addMCMove(angleMoves[3],0.1);
        relaxIntegrator.getMoveManager().addMCMove(angleMovesMulti1, 0.4);
        relaxIntegrator.getMoveManager().addMCMove(angleMovesMulti1outer, 0.4);



        //shuffleMove1.skipW=true;
//relaxIntegrator.getMoveManager().addMCMove(shuffleMove1);
        //((MCMoveStepTracker) shuffleMove1.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) angleMoves[1].getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) angleMoves[3].getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) angleMovesMulti1.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) angleMovesMulti1outer.getTracker()).setNoisyAdjustment(true);
        //((MCMoveStepTracker) shuffleMove1.getTracker()).setAcceptanceTarget(0.5);
        ((MCMoveStepTracker) angleMoves[1].getTracker()).setAcceptanceTarget(0.5);

//relaxIntegrator.getMoveManager().setFrequency(shuffleMove1,0.2);
//relaxIntegrator.getMoveManager().setFrequency(angleMoves[1], 0.8);
        ActivityIntegrate relaxActivity=new ActivityIntegrate(relaxIntegrator,10000);

        Activity activityInitRefPref = new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                angleMoves[0].skipW=true;
                angleMoves[1].skipW=true;
                for (int i = 0; i < relaxSteps*10000; i++) {

                    angleMoves[0].doTrial();
                    double chi=angleMoves[0].getChi(temperature*100);
                    if(chi<sim.getRandom().nextDouble()){
                        angleMoves[0].rejectNotify();
                    }
                    else {
                        angleMoves[0].acceptNotify();
                    }
                    angleMoves[1].doTrial();
                    chi=angleMoves[1].getChi(temperature*100);
                    if(chi<sim.getRandom().nextDouble()){
                        angleMoves[1].rejectNotify();
                    }
                    else {
                        angleMoves[1].acceptNotify();
                    }

                }
                angleMoves[0].skipW=false;
                angleMoves[1].skipW=false;
            }
        };
//sim.getController().runActivityBlocking(activityInitRefPref);
        targetCluster.setTemperature(temperature);
        IntegratorMC.dodebug=false;

        if (false){
            double size = 30;
            sim.box[0].getBoundary().setBoxSize(Vector.of(size, size, size));
            sim.box[1].getBoundary().setBoxSize(Vector.of(size, size, size));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox1.setPixelUnit(new Pixel(300.0 / size));
            //displayBox0.setShowBoundary(false);
            //displayBox1.setShowBoundary(false);
            //((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            //((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);

            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    Color[] c = new Color[]{Color.BLUE, Color.RED, Color.YELLOW, Color.GREEN};
                    return c[a.getParentGroup().getIndex()];
                }
            };
            displayBox0.setColorScheme(colorScheme);
            displayBox1.setColorScheme(colorScheme);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            ClusterAbstract myTargetCluster = targetCluster;
            Activity activityRelax = new Activity() {
                @Override
                public void runActivity(Controller.ControllerHandle handle) {
                    myTargetCluster.setTemperature(100);
                    for (int i = 0; i < 10; i++) {
                        handle.yield(sim.integrators[0]::doStep);
                        handle.yield(sim.integrators[1]::doStep);
                    }
                    myTargetCluster.setTemperature(temperature);
                }
            };
            sim.getController().addActivity(activityRelax);
            sim.initRefPref(null, 10, false);

            sim.equilibrate(null, 20, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            final DisplayTextBox stepsBox0 = new DisplayTextBox();
            stepsBox0.setLabel("Reference");
            final DisplayTextBox stepsBox1 = new DisplayTextBox();
            stepsBox1.setLabel("Target");
            JLabel jLabelPanelParentGroupSteps = new JLabel("steps");
            final JPanel panelParentGroupSteps = new JPanel(new BorderLayout());
            panelParentGroupSteps.add(jLabelPanelParentGroupSteps, Constants.CompassDirection.NORTH.toString());
            panelParentGroupSteps.add(stepsBox0.graphic(), BorderLayout.WEST);
            panelParentGroupSteps.add(stepsBox1.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroupSteps, SimulationPanel.getVertGBC());
            DataSourceCountSteps dsSteps0 = new DataSourceCountSteps(sim.integrators[0]);
            DataPumpListener pumpSteps0 = new DataPumpListener(dsSteps0, stepsBox0, 1000);
            sim.integrators[0].getEventManager().addListener(pumpSteps0);
            DataSourceCountSteps dsSteps1 = new DataSourceCountSteps(sim.integrators[1]);
            DataPumpListener pumpSteps1 = new DataPumpListener(dsSteps1, stepsBox1, 1000);
            sim.integrators[1].getEventManager().addListener(pumpSteps1);
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B" + nPoints);
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup, Constants.CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());

            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();

                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
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
        long trelax=System.nanoTime();
        sim.getController().runActivityBlocking(relaxActivity);
        trelax=System.nanoTime()-trelax;
        System.out.println("time to relax: " +trelax/ 1e9);
        angleMoves[1].skipW=false;
        angleMoves[3].skipW=false;
        angleMovesMulti1.skipW=false;
        angleMovesMulti1outer.skipW=false;
        IntegratorMC.dodebug=false;
        //angleMoves[1].skipW=false;
        //System.exit(0);

        if(false){
            IAtomList atoms = sim.box[1].getLeafList();
            long tq=System.nanoTime();
            long m=500000;
            ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[1]).setCellManager(null);
            ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[1]).skipW=true;
            for(int i=0;i<m;i++){
                sim.mcMoveRotate[1].doTrial();
                sim.mcMoveRotate[1].acceptNotify();

            }
            tq=System.nanoTime()-tq;
            System.out.println("time to rotate: " +tq/ (double) m);


        }




        if(false){
            IMoleculeList molecules = sim.box[1].getMoleculeList();
            long tq=System.nanoTime();
            long m=50000;
            for(int i=0;i<m;i++){
                pTarget.energy(molecules.get(0),molecules.get(1));
                pc1.computeAll(false);
                sim.box[1].trialNotify();
                sim.box[1].acceptNotify();
            }
            tq=System.nanoTime()-tq;
            System.out.println("time to compute energy: " +tq/ (double) m);



        }
        IntegratorMC benchmarkIntegrator= new IntegratorMC(sim.integrators[1].getPotentialCompute(), sim.getRandom(), temperature*100,sim.box[1]);

        benchmarkIntegrator.getMoveManager().addMCMove(sim.mcMoveTranslate[1]);
        benchmarkIntegrator.getMoveManager().addMCMove(sim.mcMoveRotate[1]);
        if(true){
            HistogramExpanding histogram = new HistogramExpanding(1);
            long tq=System.nanoTime();
            long m=50000;
            Vector r0=sim.box[1].getMoleculeList().get(0).getChildList().get(0).getPosition();
            Vector r1=sim.box[1].getMoleculeList().get(1).getChildList().get(0).getPosition();

            ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[1]).skipW=false;

            for(int i=0;i<m;i++){
                benchmarkIntegrator.doStep();
                double r2=r0.Mv1Squared(r1);
                histogram.addValue(Math.sqrt(r2));
            }
            tq=System.nanoTime()-tq;
            System.out.println("time to do steps: " +tq/ (double) m);
            double[]h=histogram.getHistogram();
            System.out.println("distance histogram: "+Arrays.toString(h));
            //System.exit(0);

        }
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        long initSteps = steps/40;
        if (params.fFile != null) {
            initSteps = steps;
        }
        sim.initRefPref(refFileName, initSteps);

        sim.equilibrate(refFileName, 2*initSteps);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        System.out.println("MC Move step sizes (ref)    "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[1].getStepSize())));

        sim.integratorOS.getMoveManager().setEquilibrating(false);

        ClusterAbstract tc = targetCluster;
        FileWriter fw;
        if (params.fFile != null) {
            try {
                fw = new FileWriter(params.fFile);
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            sim.integrators[1].getEventManager().addListener(new IntegratorListenerAction(new IAction() {
                @Override
                public void actionPerformed() {
                    double f = -2 * tc.value(sim.box[1]);
                    try {
                        fw.write(sim.integrators[1].getStepCount() + " " + f + "\n");
                    } catch (IOException ex) {
                        throw new RuntimeException(ex);
                    }
                }
            }));
        }
        else {
            fw = null;
        }

        sim.getController().runActivityBlocking(ai);
        long t2 = System.nanoTime();

        if (fw != null) {
            try { fw.close(); } catch (IOException ex) { throw new RuntimeException(ex); }
        }

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(refIntegral);
// Add CSV output
        String outputFile = "jack_results.csv";
        boolean fileExists = new java.io.File(outputFile).exists();

        try (FileWriter fwCsv = new FileWriter(outputFile, true)) {
            if (!fileExists) {
                fwCsv.write("numArms,armLength,coreSigma,epsCore,temperature,B2,B2_err,correlation,numSteps\n");
            }

            // Get the ratio and error from the simulation
            double[] ratioAndError = sim.dvo.getAverageAndError();
            double ratio = ratioAndError[0];
            double error = ratioAndError[1];

            // Calculate actual B2 value
            double B2 = ratio * refIntegral;
            double B2_err = error * refIntegral;

            // Get correlation from accumulator[1] (target system)
            double correlation = sim.accumulators[1].getData(sim.accumulators[1].BLOCK_CORRELATION).getValue(0);

            fwCsv.write(String.format("%d,%d,%.2f,%.2f,%.2f,%.8f,%.8f,%.8f,%d\n",
                    params.numArms,
                    params.armLength,
                    params.coreSigma,
                    params.epsCore,
                    params.temperature,
                    B2,
                    B2_err,
                    correlation,
                    params.numSteps));
            System.out.println("Results written to: " + outputFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to jack_results.csv: " + e.getMessage());
        }

        System.out.println("time: "+(t2-t1)/1e9);
    }
    // Builds polyhedral/fcc-like graft directions; normalized.


    /**
     * Parameters (defaults preserve old behavior).
     */
    public static class VirialStarParams extends ParameterBase {
        public int nPoints = 2;
        public int armLength = 4;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double refFreq = -1;
        public String fFile = null;
        public int numArms = 2;

        // Geometry & interactions
        public double coreSigma = 1.0;  // σ_core (also used in LJ)

        public double epsCore   = 1.0;  // ε_core

        // HS reference diameter options
        public boolean useCoreInHSRef = true;
        public double hsExtra = 0.0;
        public double rc=0;
    }






}

