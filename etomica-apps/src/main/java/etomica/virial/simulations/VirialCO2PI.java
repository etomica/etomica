/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.NumRootNodes;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.models.co2.P2CO2Hellmann;
import etomica.molecule.IMoleculeList;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.*;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.overlap.DataVirialOverlap.FullResult;

import javax.swing.*;
import java.awt.*;
import java.util.Map;
import java.util.Set;

//import etomica.virial.MCMoveOrientationBruteForce;

/**
 * Mayer sampling simulation for path-integral virial coefficients of CO2.
 */
public class VirialCO2PI {

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
        VirialH2PIParam params = new VirialH2PIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
//        	params.potentialLevel = levelOptions.hindePatkowski;
//        	params.subtractWhat = subOptions.iso;        	
        	params.nBeads = 4;
            params.temperature = 200;
            params.numSteps = 10000000;
            params.sigmaHSRef = 5;
            params.refFrac = -1;
//            params.doHist = true;
//            params.hackedup = true;
//            params.doHist = true;
//            params.refFrac = 0.9;
//            params.continuous = false;            
//            params.numSteps = 1000000;

        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFrac;
        
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef == -1) {
            // these correlations work fairly well over the temperature range of interest
            sigmaHSRef = 4 + 20/(10+temperatureK);
        }
        
        boolean hackedup = params.hackedup;
        if (hackedup) sigmaHSRef = 40;
        boolean continuous = params.continuous;
        if (continuous) sigmaHSRef = 7;
        int nBins = hackedup?4000:100;
        double dx = sigmaHSRef/nBins;
//        System.out.println(Mole.UNIT.fromSim(1));
//        System.exit(1);
        
        final double[] HSB = new double[8];
        final subOptions subtractWhat = params.subtractWhat;
        System.out.println("Subtract what = "+subtractWhat);
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);        
        final int beadFac = (subtractWhat == subOptions.half) ? 2 : 1;
        if (subtractWhat == subOptions.half) {
        	if (nb <= 2) throw new RuntimeException("oops!");        	
        }        
        
        if (subtractWhat == subOptions.half) {
            System.out.println("CO2 Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
            System.out.println("Calculating difference between "+nb/beadFac+" and "+nb+" beads");
        }
        else System.out.println("CO2 Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");

        if (subtractWhat == subOptions.semiClassical) {
        	System.out.println("computing difference from semiclassical");
        }
        
        final int nBeads = nb;
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
		
        final Space space = Space3D.getInstance();
        
        final double temperature = Kelvin.UNIT.toSim(temperatureK);

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        
        final P2CO2Hellmann p2c = new P2CO2Hellmann(space, P2CO2Hellmann.Parameters.B);
        
        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2c, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }

        P2CO2Hellmann.P2CO2SC p2SemiClassical = p2c.makeSemiclassical(temperature);
                
//        final MayerGeneral fTargetSemiClassical = new MayerGeneral(p2SemiClassical);

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

        boolean doFlex = nPoints > 2;
        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, doFlex);
        flexDiagrams.setDoMinimalMulti(true);
        flexDiagrams.setDoMinimalBC(true);
        flexDiagrams.setDoReeHoover(true);
        flexDiagrams.setDoShortcut(true);
        ClusterAbstract targetCluster = flexDiagrams.makeVirialCluster(fTarget, null, false);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
        final ClusterSum[] targetSubtract = new ClusterSum[(subtractWhat == subOptions.half) ? beadFac : 1];
        final ClusterSum fullTargetCluster;

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];
        int[] targetDiagramNumbers = new int[0];

        if (subtractWhat != subOptions.none) {
            fullTargetCluster = (ClusterSum)targetCluster;
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            for (int i=0; i<targetSubtract.length; i++) {
                if  (subtractWhat == subOptions.half) {
                    targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetSkip[i]});
                }                    
//                    if (subtractWhat == subOptions.semiClassical) {
//                    	targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{(fTargetSemiClassical)});
//                    }
            }

            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
            
            ClusterSumShell[] targetDiagramsPlus = flexDiagrams.makeSingleVirialClusters(fullTargetCluster, null, fTarget);
            ClusterSumShell[][] targetDiagramsMinus = new ClusterSumShell[targetDiagramsPlus.length][0];
            for (int j=0; j<targetDiagramsMinus.length; j++) {
                targetDiagramsMinus[j] = new ClusterSumShell[targetSubtract.length];
            }
            for (int i=0; i<targetSubtract.length; i++) {
                ClusterSumShell[] foo = flexDiagrams.makeSingleVirialClusters(targetSubtract[i], null, fTarget);
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
            targetDiagrams = flexDiagrams.makeSingleVirialClusters((ClusterSum)targetCluster, null, fTarget);
        }
        IsBiconnected isBi = new IsBiconnected();
        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = flexDiagrams.getMSMCGraphs(true, false);
        Map<Graph,Graph> cancelMap = flexDiagrams.getCancelMap();
        int iGraph = 0;
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ed = new DeleteEdgeParameters(flexDiagrams.eBond);
        for (Graph g : singleGraphs) {
            if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.mBond)) continue;
            if (g.nodeCount() > 3 && isBi.check(g)) {
                if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.eBond)) continue;
                System.out.print(" ("+g.coefficient()+") "+g.nodeCount()+"bc");
                targetDiagramNumbers[iGraph] = -g.nodeCount();
            }
            else {
                Graph gf = edgeDeleter.apply(g, ed);
                System.out.print(" ("+g.coefficient()+") "+gf.getStore().toNumberString());
                targetDiagramNumbers[iGraph] = Integer.parseInt(gf.getStore().toNumberString());
            }
            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                Graph gf = edgeDeleter.apply(cancelGraph, ed);
                System.out.print(" - "+gf.getStore().toNumberString());
            }
            System.out.println();
            iGraph++;
        }
        System.out.println();
        Set<Graph> disconnectedGraphs = flexDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0) {
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
        AtomTypeOriented atype = new AtomTypeOriented(new ElementSimple("CO2", Carbon.INSTANCE.getMass() + 2 * Oxygen.INSTANCE.getMass()), space);
        SpeciesSpheresHetero species = new SpeciesSpheresHetero(space, new AtomTypeOriented[]{atype}) {
            protected IAtom makeLeafAtom(AtomType leafType) {
                double bl = 2*p2c.getPos(5);
                return new AtomHydrogen(space, (AtomTypeOriented) leafType, bl);
            }
        };
        species.setChildCount(new int [] {nBeads});
        species.setConformation(new ConformationLinear(space, 0));
        
        

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints+(doFlex?1:0)}, temperature, new ClusterAbstract[]{refCluster, targetCluster},
                 targetDiagrams, new ClusterWeight[]{refSampleCluster,targetSampleCluster}, false);

//        We use ANIntragroupExchange here by purpose even though we are not doing exchange
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

//        if (hackedup) {
//            ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//            ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setTunable(false);
//            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).hackedup = true;
//            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).dx = dx;
//            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).continuous = continuous;
//        }
        
        // rotation is a bit pointless when we can regrow the chain completely
        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        
        System.out.println("regrow full ring");
        MCMoveClusterRingRegrow ring0 = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*atype.getMass()*temperature);
        ring0.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));
        
        MCMoveClusterRingRegrow ring1 = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        ring1.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));
        
        boolean fixedPosition = false;
        if (!fixedPosition) {
            System.out.println("variable position");
            sim.integrators[0].getMoveManager().addMCMove(ring0);
            sim.integrators[1].getMoveManager().addMCMove(ring1);
        }
        else {
            System.out.println("fixed position");
        }
        MCMoveClusterRingRegrowOrientation move0 = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);        
        MCMoveClusterRingRegrowOrientation move1 = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
        move0.setStiffness(temperature, Oxygen.INSTANCE.getMass());
        move1.setStiffness(temperature, Oxygen.INSTANCE.getMass());
//        System.out.println(2*nBeads*Math.PI/(lambda*lambda)+" "+2*move0.getStiffness()+" "+2*move1.getStiffness());
        boolean fixedOrientation = false;
        if (!fixedOrientation) {
            System.out.println("variable orientation");
        	sim.integrators[0].getMoveManager().addMCMove(move0);
        	sim.integrators[1].getMoveManager().addMCMove(move1);
        }        
        else {
            System.out.println("fixed orientation");
        }
//        MCMoveOrientationBruteForce orBF0 = new MCMoveOrientationBruteForce(sim.getRandom(), space, temperature);
//        MCMoveOrientationBruteForce orBF1 = new MCMoveOrientationBruteForce(sim.getRandom(), space, temperature);
//
//        ((MCMoveStepTracker)orBF0.getTracker()).setNoisyAdjustment(true);        
//        orBF0.setStepSize(0.001);
//        ((MCMoveStepTracker)orBF1.getTracker()).setNoisyAdjustment(true);
//        orBF1.setStepSize(0.001);
//        sim.integrators[0].getMoveManager().addMCMove(orBF0);
//        sim.integrators[1].getMoveManager().addMCMove(orBF1);
        


        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        if (subtractWhat != subOptions.none) {
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
            double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
            if (pi == 0) throw new RuntimeException("initialization failed");
            sim.box[1].acceptNotify();
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
            double vSize = 5;
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
            diameterManager.setDiameter(type, 0.02+1.0/nBeads);
            displayBox1.setDiameterHash(diameterManager);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
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
            refFileName += "_"+tempString+"_"+nBeads;
            if (subtractWhat == subOptions.half) refFileName += "H";
            if (subtractWhat == subOptions.none) refFileName += "N";
            if (subtractWhat == subOptions.semiClassical) refFileName += "SC";
        }
        long t1 = System.currentTimeMillis();
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);

        // this can't be done after equilibration.  ClusterSumShell needs at least
        // one accepted move before it can collect real data.  we'll reset below
        MeterVirial meterDiagrams = new MeterVirial(targetDiagrams);
        meterDiagrams.setBox(sim.box[1]);
        AccumulatorAverageCovariance accumulatorDiagrams = null;
        if (targetDiagrams.length > 0) {
            // if we have 1 diagram, we don't need covariance, but it won't actually cause problems
            accumulatorDiagrams = new AccumulatorAverageCovariance(steps);
            DataPumpListener pumpDiagrams = new DataPumpListener(meterDiagrams, accumulatorDiagrams);
            sim.integrators[1].getEventManager().addListener(pumpDiagrams);
        }

        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        
        if (accumulatorDiagrams != null) {
            accumulatorDiagrams.reset();
        }

        // make the accumulator block size equal to the # of steps performed for each overlap step.
        // make the integratorOS aggressive so that it runs either reference or target
        // then, we'll have some number of complete blocks in the accumulator
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        DataSourceScalar dsDisp = new DataSourceScalar("foo", Null.DIMENSION) {
            final Vector r = space.makeVector();
            public double getDataAsScalar() {
                IMoleculeList mols = sim.box[0].getMoleculeList();
                double sum = 0;
                int n = 0;
                for (int i = 0; i<mols.size(); i++) {
                    IAtomList atoms = mols.get(i).getChildList();
                    r.E(0);
                    for (int j = 0; j<atoms.size(); j++) {
                        r.PE(atoms.get(j).getPosition());
                    }
                    r.TE(1.0/atoms.size());
                    for (int j = 0; j<atoms.size(); j++) {
                        sum += r.Mv1Squared(atoms.get(j).getPosition());
                    }
                    n += atoms.size();
                }
                sum /= n;
                return sum;
            }
        };
        AccumulatorAverageFixed accDisp = new AccumulatorAverageFixed(10);
        DataPumpListener pumpDisp = new DataPumpListener(dsDisp, accDisp);
        sim.integrators[0].getEventManager().addListener(pumpDisp);

        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
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
        if (!isCommandline) {
            // if interactive, print intermediate results
            final double refIntegralF = refIntegral;
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.getController().getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    double a = sim.dvo.getOverlapAverage();
                    FullResult foo = sim.dvo.getFullRatioResultForAlpha(a);
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
//                    System.out.println(String.format("   alpha: %12.6e  ref: %12.6e  %8.2e   tar: %12.6e  %8.2e", a, foo.refAvg, foo.refErr, foo.targetAvg, foo.targetErr));
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
                        if ((sim.integratorOS.getStepCount()*10) % sim.getController().getMaxSteps() != 0) return;
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

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), 1000);
        long t2 = System.currentTimeMillis();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
//                    System.out.println(xValues[i]+" "+(-2*h[i]+1)+" "+Math.exp(-u/temperature));
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        System.out.println("Target Ring acceptance "+ring1.getTracker().acceptanceRatio());
        System.out.println("Target Ring orientation acceptance "+move1.getTracker().acceptanceRatio());
        System.out.println("Reference Ring orientation acceptance "+move0.getTracker().acceptanceRatio());

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
                System.out.print("diagram "+(-targetDiagramNumbers[i])+"bc ");
            }
            else {
                System.out.print("diagram "+targetDiagramNumbers[i]+" ");
            }
            // average is vi/|v| average, error is the uncertainty on that average
            // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)
            double ivar = dataCov.getValue((i+1)*nTotal+(i+1));
            double ocor = ivar*oVar == 0 ? 0 : dataCov.getValue(nTotal*(i+1)+nTotal-1)/Math.sqrt(ivar*oVar);
            System.out.print(String.format("average: %20.15e  error: %10.15e  ocor: %6.4f", dataAvg.getValue(i+1), dataErr.getValue(i+1), ocor));
            if (targetDiagrams.length > 1) {
                System.out.print("  dcor:");
                for (int j=0; j<targetDiagrams.length; j++) {
                    if (i==j) continue;
                    double jvar = dataCov.getValue((j+1)*nTotal+(j+1));
                    double dcor = ivar*jvar == 0 ? 0 : dataCov.getValue((i+1)*nTotal+(j+1))/Math.sqrt(ivar*jvar);
                    System.out.print(String.format(" %6.4f", dcor));
                }
            }
            System.out.println();
        }
//        try {
//            ((P2HydrogenPatkowski)p2).file.close();
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
        

        DataGroup dataDisp = (DataGroup)accDisp.getData();
        double avgDisp = dataDisp.getData(accDisp.AVERAGE.index).getValue(0);
        double errDisp = dataDisp.getData(accDisp.ERROR.index).getValue(0);
        double corDisp = dataDisp.getData(accDisp.BLOCK_CORRELATION.index).getValue(0);
        
        System.out.println(String.format("position stdev: %15.9e  %10.4e  %4.2f", avgDisp, errDisp, corDisp));

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
    
    enum subOptions {
        half, semiClassical, none
    }
    
    /**
     * Inner class for parameters
     */
    public static class VirialH2PIParam extends ParameterBase {
        public int nPoints = 2;
        public int nBeads = 8;
        public double temperature = 1000.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public boolean doHist = false;
        public double sigmaHSRef = 6.0; // -1 means use equation for sigmaHSRef
        
        public boolean hackedup = false;
        public boolean continuous = false;        
        public subOptions subtractWhat = subOptions.none;
    }
}
