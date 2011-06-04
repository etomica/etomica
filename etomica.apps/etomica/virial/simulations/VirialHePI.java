package etomica.virial.simulations;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IAtomType;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ANIntergroupCoupled;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.ElementChemical;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2Harmonic;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialGroup;
import etomica.space.IVectorRandom;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.Volume;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.DoubleRange;
import etomica.util.HistogramNotSoSimple;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterDifference;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumMultibody;
import etomica.virial.ClusterSumShell;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRingRegrow;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerFunction;
import etomica.virial.MayerFunctionMolecularThreeBody;
import etomica.virial.MayerFunctionNonAdditive;
import etomica.virial.MayerFunctionSphericalThreeBody;
import etomica.virial.MayerFunctionThreeBody;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirial;
import etomica.virial.PotentialGroup3PI;
import etomica.virial.PotentialGroup3PI.PotentialGroup3PISkip;
import etomica.virial.PotentialGroupPI;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Mayer sampling simulation
 */
public class VirialHePI {

    public static void main(String[] args) {

        VirialHePIParam params = new VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        final boolean pairOnly = params.nPoints == 2 || params.pairOnly;
        double refFreq = params.refFrac;
        boolean subtractHalf = params.subtractHalf;
        final double sigmaHSRef = params.sigmaHSRef;
        final int startBeadHalfs = params.startBeadHalfs;
        final int beadFac = subtractHalf ? params.beadFac : 1;
        final int finalNumBeads = params.finalNumBeads;
        final double[] HSB = new double[8];
        if (params.nBeads>-1) System.out.println("nSpheres set explicitly");
        int nb = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean doDiff = params.doDiff;
        final boolean semiClassical = params.semiClassical;
        int origNB = nb;
        if (subtractHalf) {
            int totalHalfs = (int)Math.ceil(Math.log((nb+finalNumBeads-1)/finalNumBeads)/Math.log(beadFac));
            int totalFac = (int)Math.round(Math.pow(beadFac, totalHalfs));
            int endBeads = (nb + totalFac -1) / totalFac;
            nb = endBeads * totalFac;
            origNB = nb;
            System.out.println("He Path Integral ("+origNB+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
            nb /= (int)Math.round(Math.pow(beadFac,startBeadHalfs));
            if (nb*beadFac <= finalNumBeads) {
                throw new RuntimeException("it is unnecessary to half "+(startBeadHalfs)+" times");
            }
            if (nb <= finalNumBeads) {
                if (doDiff) {
                    if (semiClassical) {
                        System.out.println("Calculating difference between "+nb+" beads and semiclassical");
                    }
                    else {
                        System.out.println("Calculating difference between "+nb+" beads and classical");
                    }
                }
                else {
                    System.out.println("Perturbing from "+nb+" beads to hard spheres");
                }
                subtractHalf = false;
            }
            else {
                System.out.println("Calculating difference between "+nb/beadFac+" and "+nb+" beads");
            }
        }
        else {
            System.out.println("He Path Integral ("+nb+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
            if (doDiff) {
                if (semiClassical) {
                    System.out.println("computing difference from semiclassical");
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
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        final Potential2SoftSpherical p2 = new P2HePCKLJS(space);

        
        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }

        final P3CPSNonAdditiveHe p3 = new P3CPSNonAdditiveHe(space);
        PotentialGroup3PI p3TargetGroup = new PotentialGroup3PI(beadFac);
        p3TargetGroup.addPotential(p3, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }

        final MayerGeneralSpherical fTargetClassical = new MayerGeneralSpherical(p2);
        P2EffectiveFeynmanHibbs p2SemiClassical = new P2EffectiveFeynmanHibbs(space, p2);
        p2SemiClassical.setMass(heMass);
        p2SemiClassical.setTemperature(temperature);
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

        boolean doFlex = nPoints > 2 && pairOnly;
        VirialDiagrams flexDiagrams = new VirialDiagrams(nPoints, true, doFlex);
        flexDiagrams.setDoMinimalMulti(true);
        ClusterAbstract targetCluster = flexDiagrams.makeVirialCluster(fTarget, null, pairOnly ? null : f3Target);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(false);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef, eRef);
        final ClusterSum[] targetSubtract = new ClusterSum[subtractHalf ? beadFac : 1];
        final ClusterSum fullTargetCluster;

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];
        int[] targetDiagramNumbers = new int[0];

        if (doDiff || subtractHalf) {
            fullTargetCluster = (ClusterSum)targetCluster;
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            for (int i=0; i<targetSubtract.length; i++) {
                if (pairOnly) {
                    targetSubtract[i] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{subtractHalf ? fTargetSkip[i] : (semiClassical ? fTargetSemiClassical : fTargetClassical)});
                }
                else {
                    targetSubtract[i] = new ClusterSumMultibody(minusBonds, wMinus, new MayerFunction[]{subtractHalf ? fTargetSkip[i] : (semiClassical ? fTargetSemiClassical : fTargetClassical)},
                        new MayerFunctionNonAdditive[]{subtractHalf ? f3TargetSkip[i] : f3TargetClassical});
                }
            }

            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
            
            ClusterSumShell[] targetDiagramsPlus = new ClusterSumShell[0];
            targetDiagramsPlus = flexDiagrams.makeSingleVirialClusters(fullTargetCluster, null, fTarget);
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
            if (pairOnly) {
                targetDiagrams = new ClusterDifference[targetDiagramsPlus.length];
                for (int j=0; j<targetDiagramsPlus.length; j++) {
                    targetDiagrams[j] = new ClusterDifference(targetDiagramsPlus[j], targetDiagramsMinus[j]);
                }
            }
        }
        else {
            if (pairOnly) {
                targetDiagrams = flexDiagrams.makeSingleVirialClusters((ClusterSum)targetCluster, null, fTarget);
            }
        }
        if (pairOnly) {
            targetDiagramNumbers = new int[targetDiagrams.length];
            System.out.println("individual clusters:");
            Set<Graph> singleGraphs = flexDiagrams.getMSMCGraphs(true, !pairOnly);
            Map<Graph,Graph> cancelMap = flexDiagrams.getCancelMap();
            int iGraph = 0;
            for (Graph g : singleGraphs) {
                if (VirialDiagrams.graphHasEdgeColor(g, flexDiagrams.mBond)) continue;
                System.out.print(" ("+g.coefficient()+") "+g.getStore().toNumberString());
                targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());
                Graph cancelGraph = cancelMap.get(g);
                if (cancelGraph != null) {
                    System.out.print(" - "+cancelGraph.getStore().toNumberString());
                }
                System.out.println();
                iGraph++;
            }
            System.out.println();
            Set<Graph> disconnectedGraphs = flexDiagrams.getExtraDisconnectedVirialGraphs();
            if (disconnectedGraphs.size() > 0) {
                System.out.println("extra clusters:");
                HashMap<Graph,Set<Graph>> splitMap = flexDiagrams.getSplitDisconnectedVirialGraphs(disconnectedGraphs);
    
                for (Graph g : disconnectedGraphs) {
                    Set<Graph> gSplit = splitMap.get(g);
                    System.out.print(g.coefficient()+" ");
                    for (Graph gs : gSplit) {
                        if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.MBond)) {
                            System.out.print(" "+gs.nodeCount()+"b");
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
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        SpeciesSpheres species = new SpeciesSpheres(space, nBeads, new AtomTypeLeaf(new ElementChemical("He", heMass, 2)), new ConformationLinear(space, 0));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints+(doFlex?1:0)}, temperature, new ClusterAbstract[]{refCluster, targetCluster},
                 new ClusterWeight[]{refSampleCluster,targetSampleCluster}, false);

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

        if (doDiff || subtractHalf || !pairOnly) {
            AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
            IVectorRandom groupTranslationVector = (IVectorRandom)translator.getTranslationVector();
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
                    IVectorMutable v = molecules.getMolecule(i).getChildList().getAtom(1).getPosition();
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
            
            
            IAtomType type = species.getLeafType();
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
                public void actionPerformed() {
                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
                
                DataDouble data = new DataDouble();
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
            refFileName = "refpref"+nPoints+"_"+tempString+"_"+nBeads;
            if (subtractHalf) {
                refFileName += "_sh";
            }
            else if (doDiff) {
                if (semiClassical) {
                    refFileName += "_sc";
                }
                else {
                    refFileName += "_c";
                }
            }
        }
        long t1 = System.currentTimeMillis();
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);

        // this can't be done after equilibration.  ClusterSumShell needs at least
        // one accepted move before it can collect real data.  we'll reset below
        MeterVirial meterDiagrams = new MeterVirial(targetDiagrams);
        meterDiagrams.setBox(sim.box[1]);
        AccumulatorAverage accumulatorDiagrams = null;
        if (targetDiagrams.length > 0) {
            if  (targetDiagrams.length > 1) {
                accumulatorDiagrams = new AccumulatorAverageCovariance(steps);
            }
            else {
                accumulatorDiagrams = new AccumulatorAverageFixed(steps);
            }
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
        sim.integratorOS.setAgressiveAdjustStepFraction(true);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IIntegratorListener histListener = new IIntegratorListener() {
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
                double v = finalTargetCluster.value(sim.box[0]);
                double r = Math.sqrt(sim.box[0].getCPairSet().getr2(0, 1));
                hist.addValue(r, v);
                if (nPoints > 2) {
                    r = Math.sqrt(sim.box[0].getCPairSet().getr2(0, 2));
                    hist.addValue(r, v);
                    r = Math.sqrt(sim.box[0].getCPairSet().getr2(1, 2));
                    hist.addValue(r, v);
                }
            }
            
            public void integratorInitialized(IIntegratorEvent e) {
            }
        };
        if (!isCommandline) {
            // if interactive, print intermediate results
            final double refIntegralF = refIntegral;
            IIntegratorListener progressReport = new IIntegratorListener() {
                public void integratorInitialized(IIntegratorEvent e) {}
                public void integratorStepStarted(IIntegratorEvent e) {}
                public void integratorStepFinished(IIntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
            if (refFreq > 0.5) {
                IIntegratorListener histReport = new IIntegratorListener() {
                    public void integratorInitialized(IIntegratorEvent e) {}
                    public void integratorStepStarted(IIntegratorEvent e) {}
                    public void integratorStepFinished(IIntegratorEvent e) {
                        if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                        double[] xValues = hist.xValues();
                        double[] h = hist.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (!Double.isNaN(h[i])) {
                                System.out.println(xValues[i]+" "+h[i]);
                            }
                        }
                    }
                };
                sim.integratorOS.getEventManager().addListener(histReport);
            }

        }
        if (refFreq > 0.5) {
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListener);
        }

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();
        
        if (refFreq > 0.5) {
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
        System.out.println("Ring acceptance "+ring0.getTracker().acceptanceRatio()+" "+ring1.getTracker().acceptanceRatio());
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+" error: "+error);
        System.out.println("abs average: "+ratio*refIntegral+" error: "+error*refIntegral);
        double conv = Mole.UNIT.toSim(1e-24);
        System.out.println("   (cm^3/mol)^"+(nPoints-1)+"  "+ratio*refIntegral*Math.pow(conv, nPoints-1)+" "+error*refIntegral*Math.pow(conv, nPoints-1));
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData();
        IData ratioData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index);
        IData ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index);
        IData averageData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.AVERAGE.index);
        IData stdevData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.ERROR.index);
        IData correlationData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_CORRELATION.index);
        IData covarianceData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_COVARIANCE.index);
        double correlationCoef = covarianceData.getValue(1)/(stdevData.getValue(0)*stdevData.getValue(1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("reference ratio average: %20.15e error:  %10.5e  cor: %9.3e\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
        System.out.print(String.format("reference average: %20.15e stdev: %9.3e error: %9.3e cor: %6.4f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("reference overlap average: %20.15e stdev: %9.3e error: %9.3e cor: %6.4f\n",
                              averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));
        
        allYourBase = (DataGroup)sim.accumulators[1].getData();
        ratioData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index);
        ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index);
        averageData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.AVERAGE.index);
        stdevData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.STANDARD_DEVIATION.index);
        errorData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.ERROR.index);
        correlationData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_CORRELATION.index);
        covarianceData = allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_COVARIANCE.index);
        correlationCoef = covarianceData.getValue(1)/(stdevData.getValue(0)*stdevData.getValue(1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("target ratio average: %20.15e  error: %10.5e  cor: %9.3e\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
        System.out.print(String.format("target average: %20.15e stdev: %9.3e error: %9.3e cor: %6.4f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("target overlap average: %20.15e stdev: %9.3e error: %9.3e cor: %6.4f\n",
                              averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));


        for (int i=0; i<targetDiagrams.length; i++) {
            DataGroup allData = (DataGroup)accumulatorDiagrams.getData();
            IData dataAvg = allData.getData(AccumulatorAverageCovariance.StatType.AVERAGE.index);
            IData dataErr = allData.getData(AccumulatorAverageCovariance.StatType.ERROR.index);
            System.out.print("diagram "+targetDiagramNumbers[i]+" ");
            System.out.print("average: "+dataAvg.getValue(i)+" error: "+dataErr.getValue(i));
            if (targetDiagrams.length > 1) {
                System.out.print(" cov:");
                IData dataCov = allData.getData(AccumulatorAverageCovariance.StatType.BLOCK_COVARIANCE.index);
                double ivar = dataCov.getValue(i*targetDiagrams.length+i);
                for (int j=0; j<targetDiagrams.length; j++) {
                    if (i==j) continue;
                    double jvar = dataCov.getValue(j*targetDiagrams.length+j);
                    System.out.print(" "+dataCov.getValue(i*targetDiagrams.length+j)/Math.sqrt(ivar*jvar));
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

    /**
     * Inner class for parameters
     */
    public static class VirialHePIParam extends ParameterBase {
        public int nPoints = 2;
        public int nBeads = 4;
        public double temperature = 100;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public boolean doDiff = true;
        public boolean semiClassical = false;
        public boolean subtractHalf = false;
        public int startBeadHalfs = 0;
        public int beadFac = 2;
        public int finalNumBeads = 2;
        public boolean pairOnly = true;
    }
}
