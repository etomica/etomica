/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConformationGeneric;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.engine.Parser;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.*;
import etomica.integrator.*;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.parser.ParserAMBER;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesAlkane;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.PotentialComputeIntramolecular;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * Mayer sampling simulation for alkanes using the AMBER force field.
  */

public class VirialAlkaneAmber {

    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagrams flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                str += " "+gs.nodeCount()+"bc";
            }
            else {
                str += " "+gs.getStore().toNumberString();
                if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
            first = false;
        }
        return str;
    }

    public static void main(String[] args) {
        VirialAlkaneAmberParam params = new VirialAlkaneAmberParam();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

        }
        final int nPoints = params.nPoints;
        //int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        String file = params.file;
        ParserAMBER.Stuff stuff= ParserAMBER.makeStuff(file, new ParserAMBER.Options());
        SpeciesManager speciesManager = stuff.speciesManager;
        ISpecies alkane = speciesManager.getSpecies(0);
        AtomType C = alkane.getAtomType(0);
        AtomType H = alkane.getAtomType(1);
        PotentialMasterBonding potentialMasterBonding = stuff.potentialMasterBonding;
        Potential2Soft [][] potential2Soft = stuff.potential2Soft;

        IntArrayList[] bonding = new IntArrayList[speciesManager.getSpecies(0).getLeafAtomCount()];
        for(int atomIndex = 0; atomIndex < bonding.length; atomIndex++){
            bonding[atomIndex] = new IntArrayList();
        }

        for(int[][][] bondType : potentialMasterBonding.getBondingInfo().bondedPairPartners[0].values()){
            for(int atomIndex = 0; atomIndex < bondType.length; atomIndex++){
                for(int[] bonds: bondType[atomIndex]){
                    int second = bonds[0] + bonds[1] - atomIndex;
                    bonding[atomIndex].add(second);
                }
            }
        }

        System.out.println("Bonding Length" + bonding.length);
        System.out.println("Bonding[0] Size" + bonding[0].size());

        double sigmaHSRef = 5;

        Space space = Space3D.getInstance();

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        PotentialMoleculePair potentialMoleculePair = new PotentialMoleculePair(space, speciesManager);
        potentialMoleculePair.setAtomPotential(C, C, potential2Soft[0][0]);
        potentialMoleculePair.setAtomPotential(C, H, potential2Soft[0][1]);
        potentialMoleculePair.setAtomPotential(H, H, potential2Soft[1][1]);

        System.out.println("nPoints"+nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);

        MayerGeneral fTarget = new MayerGeneral(potentialMoleculePair);

        boolean alkaneFlex = nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, alkaneFlex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterSum targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
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
        ClusterAbstract refCluster = new ClusterChainHS(nPoints, fRefPos);

        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        targetDiagrams = alkaneDiagrams.makeSingleVirialClusters(targetCluster, null, fTarget);
        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = alkaneDiagrams.getMSMCGraphs(true, false);
        Map<Graph,Graph> cancelMap = alkaneDiagrams.getCancelMap();
        int iGraph = 0;
        diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
            System.out.print(iGraph+" ("+g.coefficient()+") "+g.getStore().toNumberString()); // toNumberString: its corresponding number
            targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                diagramFlexCorrection[iGraph] = true;
                Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                System.out.print(" - "+getSplitGraphString(gSplit, alkaneDiagrams, false));

            }
            System.out.println();
            iGraph++;
        }
        System.out.println();
        Set<Graph> disconnectedGraphs = alkaneDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0) {
            System.out.println("extra clusters:");

            for (Graph g : disconnectedGraphs) {
                Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(g);
                System.out.println(g.coefficient()+" "+getSplitGraphString(gSplit, alkaneDiagrams, true));
            }
            System.out.println();
        }

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesManager, alkaneFlex ? (nPoints+1) : nPoints, temperature, refCluster, targetCluster);
        sim.setRandom(new RandomMersenneTwister(3000));
        sim.setExtraTargetClusters(targetDiagrams);
        sim.init();

        if (alkaneFlex) {
            int[] constraintMap = new int[nPoints + 1];
            for (int i = 0; i < nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti) sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti) sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
        }

//        ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        sim.integratorOS.setNumSubSteps(1000);

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        //sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
        //sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);

        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), space, sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
        sim.accumulators[0].setBlockSize(1);
        PotentialComputeIntramolecular potentialComputeIntramolecular0 = new PotentialComputeIntramolecular(space, sim.box[0], speciesManager, potentialMasterBonding.getBondingInfo(), potential2Soft);
        PotentialComputeIntramolecular potentialComputeIntramolecular1 = new PotentialComputeIntramolecular(space, sim.box[1], speciesManager, potentialMasterBonding.getBondingInfo(), potential2Soft);
        PotentialMasterBonding potentialMasterBonding0 = new PotentialMasterBonding(speciesManager, sim.box[0], potentialMasterBonding.getBondingInfo());
        PotentialMasterBonding potentialMasterBonding1 = new PotentialMasterBonding(speciesManager, sim.box[1], potentialMasterBonding.getBondingInfo());
        PotentialComputeAggregate potentialComputeAggregate0 = new PotentialComputeAggregate(potentialComputeIntramolecular0, potentialMasterBonding0);
        PotentialComputeAggregate potentialComputeAggregate1 = new PotentialComputeAggregate(potentialComputeIntramolecular1, potentialMasterBonding1);

        MCMoveClusterAngle mcMoveClusterAngle0 = new MCMoveClusterAngle(potentialComputeAggregate0, space, bonding, sim.getRandom(), 1);
        MCMoveClusterAngle mcMoveClusterAngle1 = new MCMoveClusterAngle(potentialComputeAggregate1, space, bonding, sim.getRandom(), 1);

        sim.integrators[0].getMoveManager().addMCMove(mcMoveClusterAngle0);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveClusterAngle1);

        MCMoveClusterTorsion mcMoveClusterTorsion0 = new MCMoveClusterTorsion(potentialComputeAggregate0, space,bonding, sim.getRandom(), 2*Math.PI);
        MCMoveClusterTorsion mcMoveClusterTorsion1 = new MCMoveClusterTorsion(potentialComputeAggregate1, space,bonding, sim.getRandom(), 2*Math.PI);

        sim.integrators[0].getMoveManager().addMCMove(mcMoveClusterTorsion0);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveClusterTorsion1);

        ((MCMoveStepTracker) mcMoveClusterAngle0.getTracker()).setNoisyAdjustment(true);
//
//        IntegratorVelocityVerletFasterer integratorMD0 = new IntegratorVelocityVerletFasterer(potentialComputeAggregate0, sim.getRandom(), 0.0005, temperature, sim.box[0]);
//        IntegratorVelocityVerletFasterer integratorMD1 = new IntegratorVelocityVerletFasterer(potentialComputeAggregate1, sim.getRandom(), 0.0005, temperature, sim.box[1]);
//
//        integratorMD0.setIsothermal(false);
//        integratorMD1.setIsothermal(false);
//        MCMoveClusterMD mcMoveClusterMD0 = new MCMoveClusterMD(integratorMD0, sim.box[0]);
//        MCMoveClusterMD mcMoveClusterMD1 = new MCMoveClusterMD(integratorMD1, sim.box[1]);
//
//        sim.integrators[0].getMoveManager().addMCMove(mcMoveClusterMD0);
//        sim.integrators[1].getMoveManager().addMCMove(mcMoveClusterMD1);

        if(false) {
            double size = (bonding.length/3 + 5) * 1.5;
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


            DiameterHashByType diameterManager = (DiameterHashByType) displayBox0.getDiameterHash();
            diameterManager.setDiameter(speciesManager.getSpecies(0).getAtomType(0), 1);
            diameterManager.setDiameter(speciesManager.getSpecies(0).getAtomType(1), 0.5);
            displayBox1.setDiameterHash(diameterManager);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);


            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
//            sim.initRefPref(null, 10, false);
//            sim.equilibrate(null, 20, false);
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
            Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints - 1});
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            averageBox.setUnit(unit);
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            errorBox.setUnit(unit);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
            return;
        }

        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);


        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));
        //System.out.println("Torsion move acceptance "+torsionMoves[0].getTracker().acceptanceRatio()+" "+
         //           torsionMoves[1].getTracker().acceptanceRatio());


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
        sim.getController().runActivityBlocking(ai);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(refIntegral);
        System.out.println("Acceptance of Torsion Move"+ mcMoveClusterTorsion0.getTracker().acceptanceProbability());

        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(sim.accumulators[1].AVERAGE.index);
        IData dataErr = allData.getData(sim.accumulators[1].ERROR.index);
        IData dataCov = allData.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        int nTotal = (targetDiagrams.length+2);
        double oVar = dataCov.getValue(nTotal*nTotal-1);

        for (int i=0; i<targetDiagrams.length; i++) {
            if (targetDiagramNumbers[i]<0) {
                System.out.print("diagram "+(-targetDiagramNumbers[i])+("bc "));
            }
            else {
                System.out.print("diagram "+targetDiagramNumbers[i]);
//                    if (fTargetDiagramNumbers[i] != 0) {
//                        System.out.print(fTargetDiagramNumbers[i]);
//                    }

                if (diagramFlexCorrection[i]) {
                    System.out.print("c");
                }
                System.out.print(" ");
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
    }

    /**
     * Inner class for parameters
     */
    public static class VirialAlkaneAmberParam extends ParameterBase {
        public int nPoints = 2;
        //public int nSpheres = 3;
        public double temperature = 298.0;// Kelvin
        public long numSteps = 1000000;
        public double refFreq = -1;
        public String file = "CCCCCC";
    }
}
