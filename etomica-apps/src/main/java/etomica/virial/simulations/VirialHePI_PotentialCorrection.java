/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.integrator.IntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ANIntergroupCoupled;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.ElementChemical;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.property.IsBiconnected;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.units.*;
import etomica.units.Dimension;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import javax.swing.*;
import java.awt.*;
import java.util.Map;
import java.util.Set;

/**
 * Mayer sampling simulation
 */
public class VirialHePI_PotentialCorrection {

    public static void main(String[] args) {
        VirialHePIParam params = new VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperature;
        long blocks = params.blocks;
        int stepsPerBlock = params.stepsPerBlock;
        long blocksEq = params.blocksEq;
        final boolean pairOnly = params.nPoints == 2 || params.pairOnly;
        double refFreq = params.refFrac;
        double sigmaHSRef = params.sigmaHSRef;
        String file = params.file;
        
        if (sigmaHSRef == -1) {
            // these correlations work fairly well over the temperature range of interest
            sigmaHSRef = 4 + 20/(10+temperatureK);
            if (!pairOnly) {
                if (nPoints == 3) {
                    sigmaHSRef -= 0.5; 
                }
                else {
                    sigmaHSRef = 4.5 + 40/(20+temperatureK);
                }
            }
        }
      
        final double[] HSB = new double[8];
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean doTotal = params.doTotal;
        if (pairOnly && doTotal) {
            throw new RuntimeException("pairOnly needs to be off to do total");
        }
       
        System.out.println("He Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
           
        System.out.println("computing difference from approximate potential");
         
        if (pairOnly) {
            System.out.println("computing pairwise contribution");
        }
        else {
            System.out.println("computing non-additive contribution");
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
        Potential2SoftSpherical p2 = new P2HePCKLJS(space);
        Potential p3 = new P3CPSNonAdditiveHe(space);
        PotentialGroupPI pTargetGroup = new PotentialGroupPI(1);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());
        PotentialGroup3PI p3TargetGroup = new PotentialGroup3PI(1);
        p3TargetGroup.addPotential(p3, new ANIntergroupCoupled(3));
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerFunctionThreeBody f3Target = new MayerFunctionMolecularThreeBody(p3TargetGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };
        
        
        Potential2SoftSpherical p2Simplified = new P2HeSimplified(space);
        P3CPSNonAdditiveHeLessSimplified p3Simplified = new P3CPSNonAdditiveHeLessSimplified(space);
        p3Simplified.setParameters(file);
    	System.out.println("simplified pair and trimer potentials used");
    	System.out.println("21 Parameters for simplified potential:");
		for (int i=0;i<p3Simplified.params.length; i++){
			System.out.println("params["+i+"] = "+p3Simplified.params[i]);
		}
        PotentialGroupPI pTargetGroupSimplified = new PotentialGroupPI(1);
        pTargetGroupSimplified.addPotential(p2Simplified, new ApiIntergroupCoupled());
        PotentialGroup3PI p3TargetGroupSimplified = new PotentialGroup3PI(1);
        p3TargetGroupSimplified.addPotential(p3Simplified, new ANIntergroupCoupled(3));
        MayerGeneral fTargetSimplified = new MayerGeneral(pTargetGroupSimplified) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerFunctionThreeBody f3TargetSimplified = new MayerFunctionMolecularThreeBody(p3TargetGroupSimplified) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads);
            }
        };

        boolean doFlex = nPoints > 2 && (pairOnly || doTotal);
        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, doFlex);
        flexDiagrams.setDoMinimalMulti(true);
        flexDiagrams.setDoMinimalBC(true);
        flexDiagrams.setDoReeHoover(true);
        flexDiagrams.setDoShortcut(true);
        ClusterAbstract targetCluster = flexDiagrams.makeVirialCluster(fTarget, pairOnly ? null : f3Target, doTotal);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];
        int[] targetDiagramNumbers = new int[0];

        
        final ClusterSum fullTargetCluster = (ClusterSum)targetCluster;
        final ClusterSum[] targetSubtract = new ClusterSum[1];
        
        ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
        double[] wMinus = fullTargetCluster.getWeights();
        
        if (pairOnly) {
            targetSubtract[0] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetSimplified});
        }
        else {
            targetSubtract[0] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{fTargetSimplified},
                new MayerFunctionNonAdditive[]{f3TargetSimplified});
        }
       

        targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
        
        ClusterSumShell[] targetDiagramsPlus = flexDiagrams.makeSingleVirialClusters(fullTargetCluster, null, fTarget);
        ClusterSumShell[][] targetDiagramsMinus  = new ClusterSumShell[targetDiagramsPlus.length][1];
        
        ClusterSumShell[] foo = flexDiagrams.makeSingleVirialClusters(targetSubtract[0], null, fTargetSimplified);
        for (int j=0; j<foo.length; j++) {
            targetDiagramsMinus[j][0] = foo[j];
        }
        
        if (pairOnly) {
            targetDiagrams = new ClusterDifference[targetDiagramsPlus.length];
            for (int j=0; j<targetDiagramsPlus.length; j++) {
                targetDiagrams[j] = new ClusterDifference(targetDiagramsPlus[j], targetDiagramsMinus[j]);
            }
        }
        

        IsBiconnected isBi = new IsBiconnected();
        if (pairOnly) {
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
                    System.out.print(g.coefficient()+" ");
                    for (Graph gs : gSplit) {
                        if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.mmBond)) {
                            System.out.print(" "+gs.nodeCount()+"m");
                        }
                        else if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                            System.out.print(" "+gs.nodeCount()+"bc");
                        }
                        else {
                            System.out.print(" "+gs.getStore().toNumberString());
                        }
                    }
                    System.out.println();
                }
                System.out.println();
            }
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
        long steps = stepsPerBlock*blocks;
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps ("+blocks+" blocks of "+stepsPerBlock+" steps)");
        System.out.println(1000+" steps per overlap-sampling block");
        SpeciesSpheres species = new SpeciesSpheres(space, nBeads, new AtomType(new ElementChemical("He", heMass, 2)), new ConformationLinear(space, 0));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints+(doFlex?1:0)}, temperature, new ClusterAbstract[]{refCluster, targetCluster},
                 targetDiagrams, new ClusterWeight[]{refSampleCluster,targetSampleCluster}, false);


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

        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

       
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
                Vector v = molecules.getMolecule(i).getChildList().getAtom(1).getPosition();
                v.TE(0.95);
            }
        }
        sim.box[1].trialNotify();
        double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        if (pi == 0) throw new RuntimeException("initialization failed");
        sim.box[1].acceptNotify();
        

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
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0/vSize));
            displayBox1.setPixelUnit(new Pixel(300.0/vSize));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);


            AtomType type = species.getLeafType();
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
//                    sim.initRefPref(null, 10);
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
            final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
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
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
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
            
            refFileName += "_potentialCorrection";
            
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
        sim.equilibrate(refFileName, blocksEq*stepsPerBlock/1000);
        
        if (accumulatorDiagrams != null) {
            accumulatorDiagrams.reset();
        }

        // make the accumulator block size equal to the # of steps performed for each overlap step.
        // make the integratorOS aggressive so that it runs either reference or target
        // then, we'll have some number of complete blocks in the accumulator
        sim.setAccumulatorBlockSize(stepsPerBlock);
        //sim.integratorOS.setNumSubSteps((int)steps);
        //sim.integratorOS.setAgressiveAdjustStepFraction(true);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        
        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IIntegratorListener histListenerRef = new IIntegratorListener() {
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
        IIntegratorListener histListenerTarget = new IIntegratorListener() {
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
            IIntegratorListener progressReport = new IIntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
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
                IIntegratorListener histReport = new IIntegratorListener() {
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
            }

        }
        if (params.doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.ai.setMaxSteps(1000);
        sim.getController().actionPerformed();
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
        IData dataAvg = allData.getData(AccumulatorAverage.AVERAGE.index);
        IData dataErr = allData.getData(AccumulatorAverage.ERROR.index);
        IData dataCov = allData.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
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
            System.out.print(String.format("average: %20.15e  error: %10.15e  ocor: %6.4f", dataAvg.getValue(i+1), dataErr.getValue(i+1), dataCov.getValue(nTotal*(i+1)+nTotal-1)/Math.sqrt(ivar*oVar)));
            if (targetDiagrams.length > 1) {
                System.out.print("  dcor:");
                for (int j=0; j<targetDiagrams.length; j++) {
                    if (i==j) continue;
                    double jvar = dataCov.getValue((j+1)*nTotal+(j+1));
                    System.out.print(String.format(" %6.4f",dataCov.getValue((i+1)*nTotal+(j+1))/Math.sqrt(ivar*jvar)));
                }
            }
            System.out.println();
        }
        
        System.out.println();

        DataGroup allData0 = (DataGroup)sim.accumulators[0].getData();
        IData dataAuto0 = allData0.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("reference autocorrelation function: "+dataAuto0.getValue(0));
        System.out.println("reference overlap autocorrelation function: "+dataAuto0.getValue(1));
        
        System.out.println();

        IData dataAuto = allData.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("target autocorrelation function: "+dataAuto.getValue(0));
        System.out.println("target overlap autocorrelation function: "+dataAuto.getValue(1));
        
        System.out.println();
        
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

    /**
     * Inner class for parameters
     */
    public static class VirialHePIParam extends ParameterBase {
        public int nPoints = 3;
        public int nBeads = 2;
        public double temperature = 2.6;   // Kelvin
        public long blocks = 1000;  //NOT overlap blocks
        public int stepsPerBlock = 1000;
        public long blocksEq=1000; //NOT overlap steps
        public double refFrac = 0.5; //not adjustment of step freqency if positive
        public boolean doHist = false;
        public double sigmaHSRef = -1; // -1 means use equation for sigmaHSRef
        public boolean pairOnly = true;
        public boolean doTotal = false;
        public String file="paramsOriginal.dat";

    }
}
