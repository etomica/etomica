/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.helium;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ANIntergroupCoupled;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Helium;
import etomica.config.ConformationLinear;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.property.NumRootNodes;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
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
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.PotentialGroup3PI.PotentialGroup3PISkip;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.simulations.SimulationVirialOverlap2;

import javax.swing.*;
import java.awt.*;
import java.util.Map;
import java.util.Set;

/**
 * Mayer sampling simulation
 */
public class VirialHePI {

    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagrams flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ed = new DeleteEdgeParameters(flexDiagrams.mmBond);
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            byte nc = gs.nodeCount();
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.mmBond)) {
                if (gs.edgeCount() < nc*(nc-1)/2) {
                    str += gs.getStore().toNumberString();
                    str += "m";
                    Graph gOnlyF = edgeDeleter.apply(edgeDeleter.apply(gs, ed), ede);
                    if (!gOnlyF.getStore().toNumberString().equals("0")) {
                        str += gOnlyF.getStore().toNumberString();
                    }
                }
                else {
                    str += " "+nc+"M";
                }
            }
            else if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
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
        VirialHePIParam params = new VirialHePIParam();
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
        final PotentialChoice pc = params.potentialChoice;
        final PotentialChoice ps = params.subtractPotential;
        final int beadFac = subtractHalf ? 2 : 1;
        final double[] HSB = new double[8];
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean doDiff = !subtractHalf && params.doDiff;
        final boolean semiClassical = params.semiClassical;
        final boolean subtractApprox = params.subtractPotential != null;
        final boolean doTotal = params.doTotal;
        if (pairOnly && doTotal) {
            throw new RuntimeException("pairOnly needs to be off to do total");
        }
        
        FlexApproach flexApproach = params.flexApproach;
        if (flexApproach != FlexApproach.FULL) {
            System.out.println("using "+flexApproach+" approach");
        }

        if (pc == PotentialChoice.APPROX) System.out.println("Calculating coefficients for approximate potential");
        else if (pc == PotentialChoice.PCKLJS) System.out.println("Calculating coefficients for PCKLJS potential");
        else if (pc == PotentialChoice.PCJS) System.out.println("Calculating coefficients for PCJS potential");
        else throw new RuntimeException("unrecognized potential");
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
                    System.out.println("computing difference from " + ps + " He");
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
        final P2HeSimplified p2Approx = new P2HeSimplified(space);
        final P2HePCKLJS p2Full = new P2HePCKLJS(space);
        final P2HePCJS p2Fuller = new P2HePCJS(space);
        final Potential2SoftSpherical p2 = (pc == PotentialChoice.APPROX) ? p2Approx : (pc == PotentialChoice.PCKLJS ? p2Full : p2Fuller);
        final Potential2SoftSpherical p2Sub = pc == null ? null : (ps == PotentialChoice.APPROX ? p2Approx : (ps == PotentialChoice.PCKLJS ? p2Full : p2Fuller));

        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }

        PotentialGroupPI pTargetSubGroup = new PotentialGroupPI(beadFac);
        pTargetSubGroup.addPotential(p2Sub, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSubSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSubSkip[i] = pTargetSubGroup.new PotentialGroupPISkip(i);
        }
        final P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);
        final P3CPSNonAdditiveHe p3Full = new P3CPSNonAdditiveHe(space);
        final IPotentialAtomicMultibody p3 = (pc == PotentialChoice.APPROX) ? p3Approx : p3Full;

        PotentialGroup3PI p3TargetGroup = new PotentialGroup3PI(beadFac);
        p3TargetGroup.addPotential(p3, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }
        final IPotentialAtomicMultibody p3Sub = (ps == PotentialChoice.APPROX ? p3Approx : p3Full);
        PotentialGroup3PI p3TargetSubGroup = new PotentialGroup3PI(beadFac);
        p3TargetSubGroup.addPotential(p3Sub, new ANIntergroupCoupled(3));

        final MayerGeneralSpherical fTargetClassical = new MayerGeneralSpherical(p2);
        Potential2Spherical p2SemiClassical = (pc == PotentialChoice.APPROX) ? p2Approx.makeQFH(temperature) : p2Full.makeQFH(temperature);
        final MayerGeneralSpherical fTargetSemiClassical = new MayerGeneralSpherical(p2SemiClassical);

        MayerGeneral[] fTargetSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            fTargetSkip[i] = new MayerGeneral(pTargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerGeneral fTargetSub = new MayerGeneral(pTargetSubGroup) {
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
        MayerFunctionThreeBody f3Target = new MayerFunctionMolecularThreeBody(p3TargetGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };
        MayerFunctionThreeBody f3TargetSub = new MayerFunctionMolecularThreeBody(p3TargetSubGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };
        boolean doFlex = (nPoints > 2 && (pairOnly || doTotal)) || nPoints > 3;
        if (flexApproach == FlexApproach.RIGID) {
            doFlex = false;
        }
        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, doFlex);
        flexDiagrams.setDoMinimalMulti(true);
        flexDiagrams.setDoMinimalBC(true);
        flexDiagrams.setDoMultiFromPair(true);
        flexDiagrams.setDoReeHoover(true);
        flexDiagrams.setDoShortcut(true);
        if (flexApproach == FlexApproach.FLEX) {
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
                            targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetSub});
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
                            targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetSub},
                                    new MayerFunctionNonAdditive[]{f3TargetSub});
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
        
        ClusterWeight targetSampleCluster = ClusterWeightAbs.makeWeightCluster(targetCluster);
        ClusterWeight refSampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);

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

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints+(doFlex?1:0)}, temperature, new ClusterAbstract[]{refCluster, targetCluster},
                 targetDiagrams, new ClusterWeight[]{refSampleCluster,targetSampleCluster}, false);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);


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

        if (false) {
            // unnecessary because our MC move regrows the chain using the
            // probability distribution appropriate for the harmonic bonds
            
            // create the intramolecular potential here, add to it and add it to
            // the potential master if needed
            PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
            // we want exp[-(pi*P/lambda^2) * sum(x^2)]
            // we set the integrator temperature=1 above, so when it does
            //   exp[-beta * U] = exp[-U]
            // so just make the spring constant whatever we need to get the above expression
            P2Harmonic p2Bond = new P2Harmonic(space, 2*Math.PI*nBeads/(lambda*lambda)*temperature);
            int[][] pairs = new int[nBeads][2];
            for (int i=0; i<nBeads-1; i++) {
                pairs[i][0] = i;
                pairs[i][1] = i+1;
            }
            pairs[nBeads-1][0] = nBeads-1;
            pairs[nBeads-1][1] = 0;
            pIntra.addPotential(p2Bond, new ApiIndexList(pairs));
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(0)});
        }

        if (false) {
            double vSize =50;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{vSize, vSize, vSize}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{vSize, vSize, vSize}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0/vSize));
            displayBox1.setPixelUnit(new Pixel(300.0/vSize));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);


            AtomType type = species.getAtomType(0);
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
            diameterManager.setDiameter(type, 1+1.0/nBeads);
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
//            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 10, false);
//                    sim.equilibrate(null, 20);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
//            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
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
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
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
        
        // if running interactively, don't use the file
        String refFileName = null;
        if (isCommandline) {
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints;
            refFileName += pairOnly ? "_2b" : "_3b";
            refFileName += "_"+tempString+"_"+nBeads;
            if (subtractHalf) {
                refFileName += "_sh";
            }
            else if (doDiff) {
                if (semiClassical) {
                    refFileName += "_sc";
                }
                else if (subtractApprox) {
                    refFileName += "_sa";
                }
                else {
                    refFileName += "_c";
                }
            }
            else {
                refFileName += "_d";
            }
            if (flexApproach == FlexApproach.RIGID) {
                refFileName += "r";
            }
            else if (flexApproach == FlexApproach.FLEX) {
                refFileName += "f";
            }
            if (pc == PotentialChoice.APPROX) {
                refFileName += "a";
            } else if (pc == PotentialChoice.PCJS) {
                refFileName += "x";
            }
        }
        long t1 = System.currentTimeMillis();
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);

        // this can't be done after equilibration.  ClusterSumShell needs at least
        // one accepted move before it can collect real data.  we'll reset below
        MeterVirial meterDiagrams = new MeterVirial(targetDiagrams);
        meterDiagrams.setBox(sim.box[1]);

        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);

        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());


        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
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
        if (!isCommandline) {
            // if interactive, print intermediate results
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
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
            if (params.doHist) {
                IntegratorListener histReport = new IntegratorListener() {
                    public void integratorInitialized(IntegratorEvent e) {}
                    public void integratorStepStarted(IntegratorEvent e) {}
                    public void integratorStepFinished(IntegratorEvent e) {
                        if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
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
            }

        }
        if (params.doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }
sim.getController().runActivityBlocking(ai);
        
        // make the accumulator block size equal to the # of steps performed for each overlap step.
        // make the integratorOS aggressive so that it runs either reference or target
        // then, we'll have some number of complete blocks in the accumulator
        long t2 = System.currentTimeMillis();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        System.out.println("Target Ring acceptance "+ring1.getTracker().acceptanceRatio());

        sim.printResults(refIntegral);

        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(sim.accumulators[1].AVERAGE.index);
        IData dataErr = allData.getData(sim.accumulators[1].ERROR.index);
        IData dataCov = allData.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        int nTotal = (targetDiagrams.length+2);
        double oVar = dataCov.getValue(nTotal*nTotal-1);
        for (int i=0; i<targetDiagrams.length; i++) {
            if (targetDiagramNumbers[i]<0) {
                System.out.print("diagram "+(-targetDiagramNumbers[i])+(pairOnly ? "bc " : "M "));
            }
            else {
                System.out.print("diagram "+targetDiagramNumbers[i]);
                if (pairOnly) {
                    if (fTargetDiagramNumbers[i] != 0) {
                        System.out.print("p"+fTargetDiagramNumbers[i]);
                    }
                }
                else {
                    System.out.print("m");
                    if (fTargetDiagramNumbers[i] != 0) {
                        System.out.print(fTargetDiagramNumbers[i]);
                    }
                }
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

        System.out.println("time: "+(t2-t1)/1000.0);
	}
    
    public static ClusterBonds[] append(ClusterBonds[] inArray, ClusterBonds[] newBonds) {
        ClusterBonds[] outArray = new ClusterBonds[inArray.length + newBonds.length];
        System.arraycopy(inArray, 0, outArray, 0, inArray.length);
        System.arraycopy(newBonds, 0, outArray, inArray.length, newBonds.length);
        return outArray;
    }

    public static double[] append(double[] inArray, double[] newWeights) {
        double[] outArray = new double[inArray.length + newWeights.length];
        System.arraycopy(inArray, 0, outArray, 0, inArray.length);
        System.arraycopy(newWeights, 0, outArray, inArray.length, newWeights.length);
        return outArray;
    }

    public enum PotentialChoice {
        APPROX,
        PCKLJS,
        PCJS
    }

    public enum FlexApproach {
        FULL,  // include all diagrams
        RIGID, // include only diagrams present in the rigid formulation
        FLEX   // include only diagrams not present in the rigid formulation (flexible correction)
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHePIParam extends ParameterBase {
        // don't change these here!!!! change them in the main method!!!
        public int nPoints = 2;
        public int nBeads = 4;
        public double temperature = 223.15;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public boolean doHist = false;
        public double sigmaHSRef = -1; // -1 means use equation for sigmaHSRef
        public boolean doDiff = true;
        public boolean semiClassical = true;
        public boolean subtractHalf = false;
        public boolean pairOnly = true;
        public boolean doTotal = false;
        public PotentialChoice potentialChoice = PotentialChoice.PCKLJS;
        public PotentialChoice subtractPotential = null;
        public FlexApproach flexApproach = FlexApproach.FULL;
    }
}
