/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.helium;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.chem.elements.Helium;
import etomica.config.ConformationLinear;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.property.NumRootNodes;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterMoleculeMulti;
import etomica.virial.mcmove.MCMoveClusterRingRegrow;
import etomica.virial.simulations.SimulationVirialOverlap2;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.util.Map;
import java.util.Set;

import static etomica.virial.simulations.helium.VirialHePI.getSplitGraphString;

/**
 * Mayer sampling simulation
 */
public class VirialHePIGraphic {

    public static void main(String[] args) {
        VirialHePI.VirialHePIParam params = new VirialHePI.VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // add any custom overrides here
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        final boolean pairOnly = params.nPoints == 2 || params.pairOnly;
        double refFreq = params.refFrac;
        boolean subtractHalf = params.subtractHalf;
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef == -1) {
            // these correlations work fairly well over the temperature range of interest
            sigmaHSRef = 2.4 + 120/(100+temperatureK);
            if (!pairOnly) {
                sigmaHSRef += 0.6;
            }
        }
        final boolean calcApprox = params.potentialChoice == VirialHePI.PotentialChoice.APPROX;
        final int beadFac = subtractHalf ? 2 : 1;
        final double[] HSB = new double[8];
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean doDiff = !subtractHalf && params.doDiff;
        final boolean semiClassical = params.semiClassical;
        final boolean subtractApprox = !calcApprox && !subtractHalf && (params.subtractPotential != null);
        final boolean doTotal = params.doTotal;
        if (pairOnly && doTotal) {
            throw new RuntimeException("pairOnly needs to be off to do total");
        }

        VirialHePI.FlexApproach flexApproach = params.flexApproach;
        if (flexApproach != VirialHePI.FlexApproach.FULL) {
            System.out.println("using "+flexApproach+" approach");
        }

        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (subtractHalf) {
            System.out.println("He Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
            System.out.println("Calculating difference between "+nb/beadFac+" and "+nb+" beads");
        }
        else {
            System.out.println("He Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
            if (doDiff) {
                if (semiClassical) {
                    System.out.println("computing difference from semiclassical");
                }
                else if (subtractApprox) {
                    System.out.println("computing difference from approximate He");
                }
                else {
                    System.out.println("computing difference from classical");
                }
            }
        }
        if (pairOnly) {
            System.out.println("computing pairwise contribution");
        }
        else {
            System.out.println("computing non-additive contribution");
        }
        if (refFreq > -1) {
            System.out.println("Setting reference step fraction to "+refFreq);
        }
        final int nBeads = nb;
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);

        Space space = Space3D.getInstance();

        double heMass = 4.002602;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        final P2HeSimplified p2Approx = new P2HeSimplified();
        final P2HePCKLJS p2Full = new P2HePCKLJS();
        final IPotential2 p2 = calcApprox ? p2Approx : p2Full;

        boolean doFlex = (nPoints > 2 && (pairOnly || doTotal)) || nPoints > 3;
        if (flexApproach == VirialHePI.FlexApproach.RIGID) {
            doFlex = false;
        }

        PotentialMoleculePairPI pTarget = new PotentialMoleculePairPI(space, p2, beadFac, nPoints+(doFlex?1:0));
        PotentialMoleculePairPI.PotentialMoleculePISkip[] pTargetSkip = new PotentialMoleculePairPI.PotentialMoleculePISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTarget.new PotentialMoleculePISkip(i);
        }
        PotentialMoleculePairPI pTargetApprox = new PotentialMoleculePairPI(space, p2, beadFac, nPoints+(doFlex?1:0));
        PotentialMoleculePairPI.PotentialMoleculePISkip[] pTargetApproxSkip = new PotentialMoleculePairPI.PotentialMoleculePISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetApproxSkip[i] = pTargetApprox.new PotentialMoleculePISkip(i);
        }

        final P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);
        final IPotential3 p3 = calcApprox ? p3Approx : new P3CPSNonAdditiveHe(space);

        PotentialMolecule3PI p3Target = new PotentialMolecule3PI(space, p3, beadFac, nPoints+(doFlex?1:0));
        PotentialMolecule3PI.PotentialMolecule3PISkip[] p3TargetSkip = new PotentialMolecule3PI.PotentialMolecule3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3Target.new PotentialMolecule3PISkip(i);
        }

        PotentialMolecule3PI p3TargetApprox = new PotentialMolecule3PI(space, p3Approx, beadFac, nPoints+(doFlex?1:0));
        PotentialMolecule3PI.PotentialMolecule3PISkip[] p3TargetApproxSkip = new PotentialMolecule3PI.PotentialMolecule3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetApproxSkip[i] = p3TargetApprox.new PotentialMolecule3PISkip(i);
        }

        final MayerGeneralSpherical fTargetClassical = new MayerGeneralSpherical(p2);
        IPotential2 p2SemiClassical = calcApprox ? p2Approx.makeQFH(temperature) : p2Full.makeQFH(temperature);
        final MayerGeneralSpherical fTargetSemiClassical = new MayerGeneralSpherical(p2SemiClassical);

        MayerGeneral[] fTargetSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            fTargetSkip[i] = new MayerGeneral(pTargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        MayerGeneral fTarget = new MayerGeneral(pTarget) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerGeneral fTargetApprox = new MayerGeneral(pTargetApprox) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };

        final MayerFunctionSphericalThreeBody f3TargetClassical = new MayerFunctionSphericalThreeBody(p3);

        MayerFunctionThreeBody[] f3TargetSkip = new MayerFunctionThreeBody[beadFac];
        for (int i=0; i<beadFac; i++) {
            f3TargetSkip[i] = new MayerFunctionMolecularThreeBody(p3TargetSkip[i]) {
                public double f(IMoleculeList pair, double[] r2, double beta) {
                    return super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        MayerFunctionThreeBody f3Target = new MayerFunctionMolecularThreeBody(p3Target) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };
        MayerFunctionThreeBody f3TargetApprox = new MayerFunctionMolecularThreeBody(p3TargetApprox) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };

        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, doFlex);
        flexDiagrams.setDoMinimalMulti(true);
        flexDiagrams.setDoMinimalBC(true);
        flexDiagrams.setDoMultiFromPair(true);
        flexDiagrams.setDoReeHoover(true);
        flexDiagrams.setDoShortcut(true);
        if (flexApproach == VirialHePI.FlexApproach.FLEX) {
            flexDiagrams.setFlexCancelOnly(true);
        }
        ClusterAbstract targetCluster = flexDiagrams.makeVirialCluster(fTarget, pairOnly ? null : f3Target, doTotal);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        if (!pairOnly) {
            rigidDiagrams.setAllPermutations(true);
        }
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
        final ClusterSum[] targetSubtract = new ClusterSum[beadFac];
        final ClusterSum fullTargetCluster;

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];

        if (doDiff || subtractHalf) {
            fullTargetCluster = (ClusterSum)targetCluster;
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            for (int i=0; i<targetSubtract.length; i++) {
                if (pairOnly) {
                    if  (subtractHalf) {
                        targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetSkip[i]});
                    }
                    else {
                        if (semiClassical) {
                            targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{(fTargetSemiClassical)});
                        }
                        else if (subtractApprox) {
                            targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetApprox});
                        }
                        else {
                            targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetClassical});
                        }
                    }
                }
                else {
                    if (subtractHalf) {
                        targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetSkip[i]},
                                new MayerFunctionNonAdditive[]{f3TargetSkip[i]});
                    }
                    else {
                        if (semiClassical) {
                            targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetSemiClassical},
                                    new MayerFunctionNonAdditive[]{f3TargetClassical});
                        }
                        else if (subtractApprox) {
                            targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetApprox},
                                    new MayerFunctionNonAdditive[]{f3TargetApprox});
                        }
                        else {
                            targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetClassical},
                                    new MayerFunctionNonAdditive[]{f3TargetClassical});
                        }
                    }
                }
            }

            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);

            ClusterSum[] targetDiagramsPlus = null;
            if (pairOnly) {
                targetDiagramsPlus = flexDiagrams.makeSingleVirialClusters(fullTargetCluster, null, fTarget);
            }
            else {
                targetDiagramsPlus = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)fullTargetCluster, fTarget, f3Target);
            }
            ClusterSum[][] targetDiagramsMinus = new ClusterSum[targetDiagramsPlus.length][0];
            for (int j=0; j<targetDiagramsMinus.length; j++) {
                targetDiagramsMinus[j] = new ClusterSum[targetSubtract.length];
            }
            for (int i=0; i<targetSubtract.length; i++) {
                ClusterSum[] foo = null;
                if (pairOnly) {
                    foo = flexDiagrams.makeSingleVirialClusters(targetSubtract[i], null, fTarget);
                }
                else {
                    foo = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)targetSubtract[i], fTarget, f3Target);
                }
                for (int j=0; j<foo.length; j++) {
                    targetDiagramsMinus[j][i] = foo[j];
                }
            }
            targetDiagrams = new ClusterDifference[targetDiagramsPlus.length];
            for (int j=0; j<targetDiagramsPlus.length; j++) {
                targetDiagrams[j] = new ClusterDifference(targetDiagramsPlus[j], targetDiagramsMinus[j]);
            }
        }
        else {
            if (pairOnly) {
                targetDiagrams = flexDiagrams.makeSingleVirialClusters((ClusterSum)targetCluster, null, fTarget);
            }
            else {
                targetDiagrams = flexDiagrams.makeSingleVirialClustersMulti((ClusterSumMultibody)targetCluster, fTarget, f3Target);
            }
        }
        int[] targetDiagramNumbers = new int[targetDiagrams.length];
        int[] fTargetDiagramNumbers = new int[targetDiagrams.length];
        boolean[] diagramFlexCorrection = new boolean[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = flexDiagrams.getMSMCGraphs(true, !pairOnly);
        Map<Graph,Graph> cancelMap = flexDiagrams.getCancelMap();
        int iGraph = 0;
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters edm = new DeleteEdgeParameters(flexDiagrams.mmBond);
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        for (Graph g : singleGraphs) {
            if (!pairOnly && NumRootNodes.value(g) > 1) continue;
            byte nc = g.nodeCount();
            if ((g.nodeCount() > 3 || !pairOnly) && g.edgeCount() == nc*(nc-1)/2) {
                if (!pairOnly) {
                    System.out.print(" ("+g.coefficient()+") "+g.nodeCount()+"M");
                    targetDiagramNumbers[iGraph] = -g.nodeCount();
                }
                else if (!VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.eBond)) {
                    System.out.print(" ("+g.coefficient()+") "+g.nodeCount()+"bc");
                    targetDiagramNumbers[iGraph] = -g.nodeCount();
                }
                else {
                    continue;
                    // skip biconnected graphs with e bonds
                }
            }
            else {
                String gnStr = g.getStore().toNumberString();
                targetDiagramNumbers[iGraph] = Integer.parseInt(gnStr);
                if (!pairOnly) {
                    Set<Graph> gSplit = flexDiagrams.getSplitDisconnectedVirialGraphs(g);
                    System.out.print(" ("+g.coefficient()+") "+getSplitGraphString(gSplit, flexDiagrams, false));
                    Graph gNoEM = edgeDeleter.apply(edgeDeleter.apply(g, ede), edm);
                    fTargetDiagramNumbers[iGraph] = Integer.parseInt(gNoEM.getStore().toNumberString());
                }
                else {
                    if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.eBond)) {
                        Graph gNoE = edgeDeleter.apply(g, ede);
                        gnStr += "p"+gNoE.getStore().toNumberString();
                        fTargetDiagramNumbers[iGraph] = Integer.parseInt(gNoE.getStore().toNumberString());
                    }
                    System.out.print(" ("+g.coefficient()+") "+gnStr);
                }
            }
            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                diagramFlexCorrection[iGraph] = true;
                String gnStr = cancelGraph.getStore().toNumberString();
                Set<Graph> gSplit = flexDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);
                if (!pairOnly) {
                    // this is actually disconnected - singlyconnected
                    if (NumRootNodes.value(cancelGraph) < NumRootNodes.value(g)) {
                        targetDiagramNumbers[iGraph] = Integer.parseInt(gnStr);
                    }
                    Graph gOnlyF = edgeDeleter.apply(cancelGraph, edm);
                    gnStr += "m";
                    fTargetDiagramNumbers[iGraph] = Integer.parseInt(gOnlyF.getStore().toNumberString());
                    if (fTargetDiagramNumbers[iGraph] > 0) {
                        gnStr += gOnlyF.getStore().toNumberString();
                    }
                    System.out.print(" - "+getSplitGraphString(gSplit, flexDiagrams, false));
                }
                else {
                    System.out.print(" - "+getSplitGraphString(gSplit, flexDiagrams, false));
                }
            }
            System.out.println();
            iGraph++;
        }
        System.out.println();
        Set<Graph> disconnectedGraphs = flexDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0 && pairOnly) {
            System.out.println("extra clusters:");

            for (Graph g : disconnectedGraphs) {
                Set<Graph> gSplit = flexDiagrams.getSplitDisconnectedVirialGraphs(g);
                if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.mmBond)) {
                    Graph cancelGraph = flexDiagrams.getCancelMap().get(g);
                    if (NumRootNodes.value(cancelGraph) < NumRootNodes.value(g)) {
                        // we have disconnected - singly connected; use the singly connected (cancelling) graph
                        gSplit = flexDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);
                    }
                }
                System.out.println(g.coefficient()+" "+getSplitGraphString(gSplit, flexDiagrams, true));
            }
            System.out.println();
        }
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        double refIntegral = HSB[nPoints];

        // the cluster's temperature determines the factor multiplied in the exponential (f=e-1)
        // we want 1/(P*kT)
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // overerr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(new AtomType(Helium.INSTANCE), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints+(doFlex?1:0)}, temperature, refCluster, targetCluster);
        sim.setExtraTargetClusters(targetDiagrams);
        sim.init();
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));


        // we'll use substeps=1000 initially (to allow for better initialization)
        // and then later switch to 1000 overlap steps
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;

        if (doFlex) {
            // fix the last molecule at the origin
            int[] constraintMap = new int[nPoints+1];
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
        }

        // rotation is a bit pointless when we can regrow the chain completely
        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);

        System.out.println("regrow full ring");
        MCMoveClusterRingRegrow ring0 = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*heMass*temperature);
        ring0.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));
        MCMoveClusterRingRegrow ring1 = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        ring1.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));

        sim.integrators[0].getMoveManager().addMCMove(ring0);
        sim.integrators[1].getMoveManager().addMCMove(ring1);

        // for flexApproach = FLEX, we need to have some non-trivial conformations
        ring1.doTrial();
        ring1.acceptNotify();
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();

        if (doDiff || subtractHalf || !pairOnly) {
            AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
            Vector groupTranslationVector = translator.getTranslationVector();
            MoleculeChildAtomAction moveMoleculeAction = new MoleculeChildAtomAction(translator);
            IMoleculeList molecules = sim.box[1].getMoleculeList();
            double r = 4;
            // put the molecules in a ring around the origin, with one atom
            // from each scaled in toward the origin
            for (int i=1; i<nPoints; i++) {
                groupTranslationVector.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
                groupTranslationVector.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
                moveMoleculeAction.actionPerformed(molecules.get(i));
                if (nBeads>1) {
                    Vector v = molecules.get(i).getChildList().get(1).getPosition();
                    v.TE(0.95);
                }
            }
            sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
            // for flexApproach = FLEX, we need may need to find some happy position
            for (int i=0; i<50 && sim.box[1].getSampleCluster().value(sim.box[1]) == 0; i++) {
                sim.mcMoveTranslate[1].doTrial();
                sim.mcMoveTranslate[1].acceptNotify();
                sim.box[1].trialNotify();
                sim.box[1].acceptNotify();
            }
            double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
            if (pi == 0) throw new RuntimeException("initialization failed");
        }
        else {
            // for flexApproach = FLEX, we need to find some happy position
            while (sim.box[1].getSampleCluster().value(sim.box[1]) == 0) {
                sim.mcMoveTranslate[1].doTrial();
                sim.mcMoveTranslate[1].acceptNotify();
                sim.box[1].trialNotify();
                sim.box[1].acceptNotify();
            }
        }

        double vSize =10;
        sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{vSize, vSize, vSize}));
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
        displayBox1.setPixelUnit(new Pixel(300.0/vSize));
        displayBox1.setShowBoundary(false);
        ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
        AtomPair pair = new AtomPair();
        for (int j=0; j<nPoints+(doFlex?1:0); j++) {
            IAtomList beads = sim.box[1].getMoleculeList().get(j).getChildList();
            for (int i=0; i<nBeads; i++) {
                pair.atom0 = beads.get(i);
                int next = i+1;
                if (next==nBeads) next=0;
                pair.atom1 = beads.get(next);
                ((DisplayBoxCanvasG3DSys)displayBox1.canvas).makeBond(pair, null);
            }
        }
        IAtomList beads = sim.box[1].getLeafList();
        for (int i=0; i<nBeads; i++) {
            pair.atom0 = beads.get(i);
            pair.atom1 = beads.get(nBeads+i);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).makeBond(pair, Color.BLUE);
        }


        AtomType type = species.getAtomType(0);
        DiameterHashByType diameterManager = (DiameterHashByType)displayBox1.getDiameterHash();
        diameterManager.setDiameter(type, 0.1+1.0/nBeads);
        displayBox1.setDiameterHash(diameterManager);
        ColorScheme colorScheme = new ColorScheme() {
            public Color getAtomColor(IAtom a) {
                return Color.RED;
            }
        };
        displayBox1.setColorScheme(colorScheme);
        simGraphic.makeAndDisplayFrame();

        final DeviceSlider tSlider = new DeviceSlider(sim.getController());
        tSlider.setShowBorder(true);
        tSlider.setBorderAlignment(TitledBorder.CENTER);
        tSlider.setShowValues(true);
        tSlider.setMinimum(0);
        tSlider.setPrecision(1);
        tSlider.setMaximum(4);
        tSlider.setNMajor(4);
        tSlider.setEditValues(true);

        final ModifierTPI tmod = new ModifierTPI(heMass, nBeads, ring0, ring1, targetCluster);
        tmod.setValue(2);
        tSlider.setModifier(tmod);
        tSlider.setLabel("log10(T)");

        final DisplayTextBox tDisplay = new DisplayTextBox();
        DataDouble.DataInfoDouble tInfo = new DataInfoDouble("T(K)", Null.DIMENSION);
        tDisplay.putDataInfo(tInfo);

        final IAction updateT = new IAction() {
            DataDouble data = new DataDouble();

            public void actionPerformed() {
                data.x = Math.pow(10, tSlider.getValue());
                tDisplay.putData(data);
            }
        };
        updateT.actionPerformed();

        tSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                updateT.actionPerformed();
            }
        });

        simGraphic.add(tSlider);
        simGraphic.getPanel().controlPanel.add(tDisplay.graphic(), SimulationPanel.getVertGBC());

        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(1000);

        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }

        final DisplayTextBox averageBox = new DisplayTextBox();
        averageBox.setLabel("Average");
        final DisplayTextBox errorBox = new DisplayTextBox();
        errorBox.setLabel("Error");
        JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
        final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
        panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
        panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
        panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);

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
        IDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new etomica.units.dimensions.Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
        Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints-1});
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

    public static class ModifierTPI implements Modifier {

        protected final MCMoveClusterRingRegrow ring0, ring1;
        protected final ClusterAbstract targetCluster;
        protected final double heMass;
        protected final int nBeads;
        protected double oldValue = 0;

        public ModifierTPI(double heMass, int nBeads, MCMoveClusterRingRegrow ring0, MCMoveClusterRingRegrow ring1, ClusterAbstract targetCluster) {
            this.heMass = heMass;
            this.nBeads = nBeads;
            this.ring0 = ring0;
            this.ring1 = ring1;
            this.targetCluster = targetCluster;
        }

        public void setValue(double newValue) {
            double temperature = Kelvin.UNIT.toSim(Math.pow(10, newValue));
            double lambda = Constants.PLANCK_H / Math.sqrt(2 * Math.PI * heMass * temperature);
            ring0.setEnergyFactor(nBeads * Math.PI / (lambda * lambda));
            ring1.setEnergyFactor(nBeads * Math.PI / (lambda * lambda));
            targetCluster.setTemperature(temperature);
            oldValue = newValue;
        }

        public double getValue() {
            // cheat.
            return oldValue;
        }

        public Dimension getDimension() {
            return Null.DIMENSION;
        }

        public String getLabel() {
            return "";
        }
    }
}
