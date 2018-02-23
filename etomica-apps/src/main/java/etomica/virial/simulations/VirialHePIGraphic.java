/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.*;
import etomica.atom.iterator.ANIntergroupCoupled;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.ElementChemical;
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
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheres;
import etomica.units.*;
import etomica.units.dimensions.*;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.PotentialGroup3PI.PotentialGroup3PISkip;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.util.Map;
import java.util.Set;

import static etomica.virial.simulations.VirialHePI.getSplitGraphString;

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
        final boolean calcApprox = params.calcApprox;
        final int beadFac = subtractHalf ? 2 : 1;
        final double[] HSB = new double[8];
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean doDiff = !subtractHalf && params.doDiff;
        final boolean semiClassical = params.semiClassical;
        final boolean subtractApprox = !calcApprox && !subtractHalf && params.subtractApprox;
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
        final P2HeSimplified p2Approx = new P2HeSimplified(space);
        final P2HePCKLJS p2Full = new P2HePCKLJS(space);
        final Potential2SoftSpherical p2 = calcApprox ? p2Approx : p2Full;

        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }

        PotentialGroupPI pTargetApproxGroup = new PotentialGroupPI(beadFac);
        pTargetApproxGroup.addPotential(p2Approx, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetApproxSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetApproxSkip[i] = pTargetApproxGroup.new PotentialGroupPISkip(i);
        }
        final P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);
        final IPotentialAtomicMultibody p3 = calcApprox ? p3Approx : new P3CPSNonAdditiveHe(space);

        PotentialGroup3PI p3TargetGroup = new PotentialGroup3PI(beadFac);
        p3TargetGroup.addPotential(p3, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }
        PotentialGroup3PI p3TargetApproxGroup = new PotentialGroup3PI(beadFac);
        p3TargetApproxGroup.addPotential(p3Approx, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetApproxSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetApproxSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }

        final MayerGeneralSpherical fTargetClassical = new MayerGeneralSpherical(p2);
        Potential2Spherical p2SemiClassical = calcApprox ? p2Approx.makeQFH(temperature) : p2Full.makeQFH(temperature);
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
        MayerGeneral fTargetApprox = new MayerGeneral(pTargetApproxGroup) {
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
        MayerFunctionThreeBody f3TargetApprox = new MayerFunctionMolecularThreeBody(p3TargetApproxGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };
        boolean doFlex = (nPoints > 2 && (pairOnly || doTotal)) || nPoints > 3;
        if (flexApproach == VirialHePI.FlexApproach.RIGID) {
            doFlex = false;
        }
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
        
        ClusterWeight targetSampleCluster = ClusterWeightAbs.makeWeightCluster(targetCluster);
        ClusterWeight refSampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // overerr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        SpeciesSpheres species = new SpeciesSpheres(space, nBeads, new AtomType(new ElementChemical("He", heMass, 2)), new ConformationLinear(space, 0));

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
                moveMoleculeAction.actionPerformed(molecules.getMolecule(i));
                if (nBeads>1) {
                    Vector v = molecules.getMolecule(i).getChildList().get(1).getPosition();
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

        double vSize =10;
        sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
        sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
        DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
        displayBox0.setPixelUnit(new Pixel(300.0/vSize));
        displayBox1.setPixelUnit(new Pixel(300.0/vSize));
        displayBox0.setShowBoundary(false);
        displayBox1.setShowBoundary(false);
        ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
        ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
        AtomPair pair = new AtomPair();
        for (int j=0; j<nPoints+(doFlex?1:0); j++) {
            IAtomList beads = sim.box[1].getMoleculeList().getMolecule(j).getChildList();
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


        AtomType type = species.getLeafType();
        DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
        diameterManager.setDiameter(type, 0.1+1.0/nBeads);
        displayBox1.setDiameterHash(diameterManager);
        ColorScheme colorScheme = new ColorScheme() {
            public Color getAtomColor(IAtom a) {
                return Color.RED;
            }
        };
        displayBox0.setColorScheme(colorScheme);
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
}
