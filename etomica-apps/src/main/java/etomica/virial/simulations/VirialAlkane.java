/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.util.GraphNumber;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesAlkane;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.ClusterChainHS;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.ClusterSumShell;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.mcmove.MCMoveClusterMoleculeHSChain;
import etomica.virial.mcmove.MCMoveClusterMoleculeMulti;
import etomica.virial.mcmove.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.mcmove.MCMoveClusterTorsionMulti;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Mayer sampling simulation for alkanes using the TraPPE-United Atoms force field.
 *
 *   M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase
 *   Equilibria. 1. United-Atom Description of n-Alkanes," J. Phys. Chem. B
 *   102, 2569-2577 (1998)
 *   
 *   
 *   modified from VirialAlkaneFlex2 so that eovererr can be used
 *   March 2013
 */


public class VirialAlkane {

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
        VirialSiepmannSpheresParam params = new VirialSiepmannSpheresParam();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 3; // B order
            params.nSpheres = 3; // length of alkane
            params.temperature = 300;
            params.numSteps = 100000000;
        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double sigmaCH2 = 3.95;
        double sigmaCH3 = nSpheres > 1 ? 3.75 : 3.73;
        double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;

        Space space = Space3D.getInstance();

        SpeciesGeneral species = SpeciesAlkane.create(nSpheres);
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(nSpheres+"-mer(TraPPE-UA) B"+nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(nSpheres > 1 ? 98.0 : 148);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        boolean alkaneFlex = nSpheres > 2 && nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, alkaneFlex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterSum targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        if (nSpheres > 2) {
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
        }

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);

        if (nSpheres > 2) {
            P3BondAngle p3 = new P3BondAngle(Math.PI*114.0/180.0, Kelvin.UNIT.toSim(62500));
            List<int[]> triplets = new ArrayList<>();
            for (int i=0; i<nSpheres-2; i++) {
                triplets.add(new int[]{i,i+1,i+2});
            }
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }

        P4BondTorsion p4 = null;
        if (nSpheres > 3) {
            p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            List<int[]> quads = new ArrayList<>();
            for (int i=0; i<nSpheres-3; i++) {
                quads.add(new int[]{i,i+1,i+2,i+3});
            }
            bondingInfo.setBondingPotentialQuad(species, p4, quads);
        }

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm,
                new int[]{alkaneFlex ? (nPoints+1) : nPoints},temperature, refCluster, targetCluster);
        sim.setExtraTargetClusters(targetDiagrams);
        sim.setDoWiggle(nSpheres > 2);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.setRandom(new RandomMersenneTwister(2));
        sim.init();

//        sim.integrators[0].getMoveManager().remove
//        MCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
        sim.accumulators[0].setBlockSize(1);

        if (alkaneFlex) {
            int[] constraintMap = new int[nPoints+1];
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
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
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        AtomType typeCH3 = species.getAtomType(0);
        pTarget.setAtomPotential(typeCH3, typeCH3, p2CH3);
        AtomType typeCH2 = nSpheres > 2 ?  species.getAtomType(1) : null;
        if (nSpheres>2) {
            pTarget.setAtomPotential(typeCH2, typeCH2, p2CH2);
            pTarget.setAtomPotential(typeCH2, typeCH3, p2CH2CH3);
        }

        sim.integratorOS.setNumSubSteps(1000);

        // create the intramolecular potential here, add to it and add it to
        // the potential master if needed
        MCMoveClusterTorsionMulti[] torsionMoves = null;
        
        if (nSpheres > 3) {
            torsionMoves = new MCMoveClusterTorsionMulti[2];
            torsionMoves[0] = new MCMoveClusterTorsionMulti(sim.integrators[0].getPotentialCompute(), space, sim.getRandom(), p4, 40);
            torsionMoves[0].setBox(sim.box[0]);
            torsionMoves[0].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMoves[0]);
            torsionMoves[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialCompute(), space, sim.getRandom(), p4, 40);
            torsionMoves[1].setBox(sim.box[1]);
            torsionMoves[1].setTemperature(temperature);
            sim.integrators[1].getMoveManager().addMCMove(torsionMoves[1]);
        }

        if(false) {
            double size = (nSpheres + 5) * 1.5;
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
            System.out.println("Flag");

            DiameterHashByType diameterManager = (DiameterHashByType) displayBox0.getDiameterHash();
            if (nSpheres>2) diameterManager.setDiameter(typeCH2, 0.2 * sigmaCH2);
            diameterManager.setDiameter(typeCH3, 0.2 * sigmaCH3);
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
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20, false);
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

        long t1 = System.nanoTime();
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
        System.out.println("MC Move step sizes (ref)    "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[0].getStepSize()))+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[1].getStepSize()))+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(ai);
        long t2 = System.nanoTime();

        if (nSpheres > 3) {
            System.out.println("Torsion move acceptance "+torsionMoves[0].getTracker().acceptanceRatio()+" "+
                    torsionMoves[1].getTracker().acceptanceRatio());
        }

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        String[] extraNames = new String[targetDiagrams.length];
        for (int i=0; i<targetDiagrams.length; i++) {
            String n = "";
            if (targetDiagramNumbers[i] < 0) {
                n = "diagram " + (-targetDiagramNumbers[i]) + "bc";
            } else {
                n = "diagram " + targetDiagramNumbers[i];

                if (diagramFlexCorrection[i]) {
                    n += "c";
                }
            }
            extraNames[i] = n;
        }
        sim.printResults(refIntegral, extraNames);
        System.out.println("time: "+(t2-t1)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialSiepmannSpheresParam extends ParameterBase {
        public int nPoints = 4;
        public int nSpheres = 3;
        public double temperature = 298.0;// Kelvin
        public long numSteps = 1000000;
        public double refFreq = -1;
    }
}
