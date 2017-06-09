package etomica.virial.simulations;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.data.*;
import etomica.integrator.IntegratorListener;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.math.function.IFunction;
import etomica.atom.IMoleculeList;
import etomica.api.IPotential;
import etomica.space.Vector;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.chem.elements.ElementSimple;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.potential.P2SpheroPolyhedron;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesPolyhedron;
import etomica.units.Null;
import etomica.math.DoubleRange;
import etomica.data.histogram.HistogramReweightedData;
import etomica.data.histogram.HistogramSimple;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.CalcFFT;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterChainHS;
import etomica.virial.ClusterSinglyConnected;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ClusterWeightUmbrella;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.ClusterWheatleyPartitionScreening;
import etomica.virial.MCMoveClusterAtomHSChain;
import etomica.virial.MCMoveClusterAtomHSRing;
import etomica.virial.MCMoveClusterAtomHSTree;
import etomica.virial.MCMoveClusterAtomQ;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralAtomic;
import etomica.virial.MeterVirialBD;
import etomica.virial.cluster.Standard;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialPolyhedra {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.shapes = new String[]{"Cube"};
            params.nShapes = new int[]{2};
            params.nPoints = 2;
            params.numSteps = 1000000L;
            params.ref = VirialHSParam.RING_TREE;
            params.doHist = false;
            params.chainFrac = 1;
            params.ringFrac = 0.7;
            params.doResize = false;
        }

        final String[] shapes = params.shapes;
        final int[] nShapes = params.nShapes;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;
        boolean doHist = params.doHist;
        final double chainFrac = params.chainFrac;
        final double ringFrac = params.ringFrac;
//        final int nSpheres = params.nSpheres;
//        final int nTO = params.nTO;
//        final int nCubes = params.nCubes;
        final boolean doResize = params.doResize;

        System.out.println("Polyhedra singly-connected sampling B"+nPoints);
//        if (nSpheres>0) {
//            System.out.println(nSpheres+" spheres");
//        }
//        if (nTO>0) {
//            System.out.println(nTO+" truncated octahedrons");
//        }
//        if (nCubes>0) {
//            System.out.println(nCubes+" cubes");
//        }
        for (int i=0; i<shapes.length; i++) {
            System.out.println(nShapes[i]+" "+shapes[i]);
        }
        

        Space space = Space3D.getInstance();

        P2SpheroPolyhedron p2 = new P2SpheroPolyhedron(space);
        MayerFunction fTarget = new MayerGeneralAtomic(p2);

        ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fTarget);
//        double v0 = 32;
//        double vs = Math.PI/6;
//        double lr = doResize ? Math.pow(vs/v0, 1.0/3.0) : 1;
//
//        List<Vector> verticesTO = new ArrayList<Vector>(24);
//        for (int a=-1; a<=1; a+=2) {
//            for (int b=-2; b<=2; b+=4) {
//                IVectorMutable v = space.makeVector();
//                v.E(new double[]{0,a*lr,b*lr});
//                verticesTO.add(v);
//                v = space.makeVector();
//                v.E(new double[]{0,b*lr,a*lr});
//                verticesTO.add(v);
//                v = space.makeVector();
//                v.E(new double[]{a*lr,0,b*lr});
//                verticesTO.add(v);
//                v = space.makeVector();
//                v.E(new double[]{b*lr,0,a*lr});
//                verticesTO.add(v);
//                v = space.makeVector();
//                v.E(new double[]{a*lr,b*lr,0});
//                verticesTO.add(v);
//                v = space.makeVector();
//                v.E(new double[]{b*lr,a*lr,0});
//                verticesTO.add(v);
//            }
//        }
//        List<Vector> verticesSphere = new ArrayList<Vector>(0);
//        List<Vector> verticesCube = new ArrayList<Vector>(8);
        
//        v0 = 8;
//        lr = doResize ? Math.pow(vs/v0, 1.0/3.0) : 1;
//        for (int a=-1; a<=1; a+=2) {
//            for (int b=-1; b<=1; b+=2) {
//                for (int c=-1; c<=1; c+=2) {
//                    IVectorMutable v = space.makeVector();
//                    v.E(new double[]{a*lr,b*lr,c*lr});
//                    verticesCube.add(v);
//                }
//            }
//        }
//
//
//        Species speciesTO = new SpeciesPolyhedron(space, verticesTO, 0.0, new ElementSimple("TO"));
//        Species speciesCube = new SpeciesPolyhedron(space, verticesCube, 0.0, new ElementSimple("C"));

//        final double sigmaTO = 2*((AtomTypeSpheroPolyhedron)speciesTO.getAtomType(0)).getOuterRadius();
//        double sigmaTOInner = 2*((AtomTypeSpheroPolyhedron)speciesTO.getAtomType(0)).getInnerRadius();
//        final double sigmaCube = 2*((AtomTypeSpheroPolyhedron)speciesCube.getAtomType(0)).getOuterRadius();
//        final double sigmaHS = (doResize || nSpheres>0) ? 1.0 : (nTO>0 ? sigmaTO : sigmaCube);
//        if (nTO>0) System.out.println("TO inner: "+sigmaTOInner+"  outer: "+sigmaTO);
//        if (nCubes>0) System.out.println("cube size: "+lr*2);
//        if (nSpheres>0) System.out.println("sphere sigma: "+sigmaHS);
//        Species speciesSphere = new SpeciesPolyhedron(space, verticesSphere, sigmaHS/2, new ElementSimple("S"));

//        double shsref = 1;
//        if (nCubes>0 && sigmaCube>shsref) shsref = sigmaCube;
//        if (nTO>0 && sigmaTO>shsref) shsref = sigmaTO;
//        final double sigmaHSRef = shsref;
        
        SpeciesPolyhedron[] allSpecies = new SpeciesPolyhedron[shapes.length];
        double shsref = 0;
        for (int i=0; i<shapes.length; i++) {
            List<Vector> vertices = ShapeParser.doParse("shape/"+shapes[i]+".dat", space).vertices;
            allSpecies[i] = new SpeciesPolyhedron(space, vertices, 0.0, new ElementSimple("P"+i));
            double s = 2*((AtomTypeSpheroPolyhedron)allSpecies[i].getAtomType(0)).getOuterRadius();
            shsref = shsref < s ? s : shsref;
        }
        final double sigmaHSRef = shsref;
        System.out.println("reference diameter: "+sigmaHSRef);

        double vhs = (4.0/3.0)*Math.PI*sigmaHSRef*sigmaHSRef*sigmaHSRef;
        MayerFunction fRefPos = new MayerFunction() {

            public void setBox(Box box) {}
            public IPotential getPotential() {return null;}

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef*sigmaHSRef ? 1 : 0;
            }
        };

        targetCluster.setTemperature(1.0);

        ClusterAbstract refCluster = null;
        long numDiagrams = 0;

        double ri = 0;
        if (ref == VirialHSParam.TREE) {
            System.out.println("using a tree reference");
            refCluster = new ClusterSinglyConnected(nPoints, fRefPos);
            numDiagrams = ((ClusterSinglyConnected)refCluster).numDiagrams();
            ri = numDiagrams*Math.pow(vhs, nPoints-1);
        }
        else if (ref == VirialHSParam.CHAINS) {
            System.out.println("using a chain reference");
            refCluster = new ClusterChainHS(nPoints, fRefPos);
            numDiagrams = ((ClusterChainHS)refCluster).numDiagrams();
            ri = numDiagrams*Math.pow(vhs, nPoints-1);
        }
        else if (ref == VirialHSParam.RINGS) {
            System.out.println("using a ring reference");
            refCluster = new ClusterChainHS(nPoints, fRefPos, true);
            numDiagrams = ((ClusterChainHS)refCluster).numDiagrams();
            final double dr = 0.00001;
            CalcFFT myFFT = new CalcFFT(new IFunction() {
                public double f(double x) {
                    if (Math.abs(x-sigmaHSRef) < 0.1*dr) {
                        return 0.5;
                    }
                    return x<sigmaHSRef ? 1 : 0;
                }
            }, dr*sigmaHSRef, 20);
            List<Object> strands = new ArrayList<Object>();
            strands.add(2);
            List<Integer> list1 = new ArrayList<Integer>();
            for (int i=1; i<nPoints; i++) {
                list1.add(2);
            }
            strands.add(list1);
            List<Object> oneMore = new ArrayList<Object>();
            oneMore.add(strands);
            ri = numDiagrams*myFFT.value(oneMore, true)[0][0];
        }
        else if (ref == VirialHSParam.CHAIN_TREE) {
            System.out.println("using a chain/tree reference ("+chainFrac+" chains)");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }
            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{chainFrac/(SpecialFunctions.factorial(nPoints)/2),(1-chainFrac)/numTreeDiagrams});
            ri = Math.pow(vhs, nPoints-1);
        }
        else if (ref == VirialHSParam.RING_TREE) {
            System.out.println("using a ring/tree reference ("+ringFrac+" rings)");
            ClusterChainHS cr = new ClusterChainHS(nPoints, fRefPos, true);
            long numRingDiagrams = cr.numDiagrams();
            double ringIntegral = numRingDiagrams*Standard.ringHS(nPoints)*Math.pow(sigmaHSRef, 3*(nPoints-1));
            double chainIntegral = (SpecialFunctions.factorial(nPoints)/2)*Math.pow(vhs, nPoints-1);
            ClusterChainHS crc = new ClusterChainHS(nPoints, fRefPos, 0, ringFrac/ringIntegral);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{crc, ct});
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }

            System.out.println("ring integral: "+ringIntegral);
            double treeIntegral = numTreeDiagrams*Math.pow(vhs, nPoints-1);

            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{1,(1-ringFrac)/treeIntegral});
            ri = 1;
        }
        else if (ref == VirialHSParam.RING_CHAIN_TREES) {
            System.out.println("using a ring/chain/tree reference ("+ringFrac+" rings, "+chainFrac+" chains)");
            ClusterChainHS cr = new ClusterChainHS(nPoints, fRefPos, true);
            long numRingDiagrams = cr.numDiagrams();
            double ringIntegral = numRingDiagrams*Standard.ringHS(nPoints);
            double chainIntegral = (SpecialFunctions.factorial(nPoints)/2)*Math.pow(vhs, nPoints-1);
            ClusterChainHS crc = new ClusterChainHS(nPoints, fRefPos, chainFrac/chainIntegral, ringFrac/ringIntegral);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{crc, ct});
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }

            double treeIntegral = numTreeDiagrams*Math.pow(vhs, nPoints-1);

            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{1,(1-ringFrac-chainFrac)/treeIntegral});
            ri = 1;
        }


        // (nPoints-1)! is simply not included by ClusterWheatley, so do that here.
        final double refIntegral = ri/SpecialFunctions.factorial(nPoints);
        System.out.println("reference integral: "+refIntegral);
        refCluster.setTemperature(1.0);



        System.out.println(steps+" steps");

        final SimulationVirial sim = new SimulationVirial(space, allSpecies, nShapes, 1.0,ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, new ClusterAbstract[]{targetCluster});
        sim.init();
        MeterVirialBD meter = new MeterVirialBD(sim.allValueClusters);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        AccumulatorAverageFixed accumulator = sim.accumulator;
        accumulator = new AccumulatorAverageFixed(1000);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(100000000);

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        if (ref == VirialHSParam.TREE) {
            MCMoveClusterAtomHSTree mcMoveTree = new MCMoveClusterAtomHSTree(sim.getRandom(), space, sigmaHSRef);
            sim.integrator.getMoveManager().addMCMove(new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveTree));
        }
        else if (ref == VirialHSParam.CHAINS) {
            MCMoveClusterAtomHSChain mcMoveHS = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
            sim.integrator.getMoveManager().addMCMove(new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHS));
        }
        else if (ref == VirialHSParam.RINGS) {
            MCMoveClusterAtomHSRing mcMoveHS = new MCMoveClusterAtomHSRing(sim.getRandom(), space, sigmaHSRef);
            sim.integrator.getMoveManager().addMCMove(new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHS));
        }
        else if (ref == VirialHSParam.CHAIN_TREE) {
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), space, sigmaHSRef);
            MCMoveClusterAtomQ mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHST);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, 1-chainFrac);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
            mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHSC);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, chainFrac);
        }
        else if (ref == VirialHSParam.RING_TREE) {
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), space, sigmaHSRef);
            MCMoveClusterAtomQ mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHST);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, 1-ringFrac);
            MCMoveClusterAtomHSRing mcMoveHSCR = new MCMoveClusterAtomHSRing(sim.getRandom(), space, sigmaHSRef);
            mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHSCR);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, ringFrac);
        }
        else if (ref == VirialHSParam.RING_CHAIN_TREES) {
            MCMoveClusterAtomHSRing mcMoveHSR = new MCMoveClusterAtomHSRing(sim.getRandom(), space, sigmaHSRef);
            MCMoveClusterAtomQ mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHSR);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, ringFrac);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
            mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHSC);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, chainFrac);
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), space, sigmaHSRef);
            mcmcaq = new MCMoveClusterAtomQ(sim.getRandom(), space, mcMoveHST);
            sim.integrator.getMoveManager().addMCMove(mcmcaq);
            sim.integrator.getMoveManager().setFrequency(mcmcaq, 1-ringFrac-chainFrac);
        }

        final HistogramReweightedData histTarg = new HistogramReweightedData(100, new DoubleRange(0, nPoints/2.0));
        final AccumulatorHistogram accHistTarg = new AccumulatorHistogram(histTarg);
        accHistTarg.setPushInterval(1000);
        final HistogramSimple histRef = new HistogramSimple(100, new DoubleRange(0, nPoints/2.0));
        final AccumulatorHistogram accHistRef = new AccumulatorHistogram(histRef);
        accHistRef.setPushInterval(1000);
        final HistogramSimple histRef0 = new HistogramSimple(100, new DoubleRange(0, nPoints/2.0));
        final AccumulatorHistogram accHistRef0 = new AccumulatorHistogram(histRef0);
        accHistRef0.setPushInterval(1000);
        final HistogramReweightedData histTargRingy = new HistogramReweightedData(100, new DoubleRange(0, nPoints/2.0));
        final AccumulatorHistogram accHistTargRingy = new AccumulatorHistogram(histTargRingy);
        accHistTargRingy.setPushInterval(1000);
        final HistogramSimple histRefRingy = new HistogramSimple(100, new DoubleRange(0, nPoints/2.0));
        final AccumulatorHistogram accHistRefRingy = new AccumulatorHistogram(histRefRingy);
        accHistRefRingy.setPushInterval(1000);
        DataInfoDoubleArray diTarg = new DataInfoDoubleArray("hist", Null.DIMENSION, new int[]{2});
        DataInfoDouble diRef = new DataInfoDouble("hist", Null.DIMENSION);
        DataTag ditag = new DataTag();
        diTarg.addTag(ditag);
        diRef.addTag(ditag);
        accHistTarg.putDataInfo(diTarg);
        accHistRef.putDataInfo(diRef);
        accHistRef0.putDataInfo(diRef);
        accHistTargRingy.putDataInfo(diTarg);
        accHistRefRingy.putDataInfo(diRef);
        final ClusterAbstract finalTargetCluster = targetCluster;
        final ClusterAbstract finalRefCluster = refCluster;
        IntegratorListener histListener = new IntegratorListener() {
            DataDoubleArray dataTarg = new DataDoubleArray(2);
            DataDouble dataRef = new DataDouble();
            Vector com = Space3D.getInstance().makeVector();
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double v = Math.abs(finalTargetCluster.value(sim.box));
//                if (v == 0) return;
                v /= finalRefCluster.value(sim.box);
                double[] x = dataTarg.getData();
                x[1] = v;
                com.E(0);
                for (int i=0; i<nPoints; i++) {
                    com.PE(sim.box.getLeafList().getAtom(i).getPosition());
                }
                com.TE(1.0/nPoints);
                for (int i=0; i<nPoints; i++) {
                    double r2 = sim.box.getLeafList().getAtom(i).getPosition().Mv1Squared(com);
                    x[0] = Math.sqrt(r2);
                    if (v!= 0) accHistTarg.putData(dataTarg);
                    dataRef.x = x[0];
                    accHistRef.putData(dataRef);
                }
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        if (doHist) sim.integrator.getEventManager().addListener(histListener);
        IntegratorListener histListenerRingy = new IntegratorListener() {
            DataDoubleArray dataTarg = new DataDoubleArray(2);
            DataDouble dataRef = new DataDouble();
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double v = Math.abs(finalTargetCluster.value(sim.box));
//                if (v == 0) return;
                v /= finalRefCluster.value(sim.box);
                double[] x = dataTarg.getData();
                x[1] = v;
                for (int i=0; i<nPoints; i++) {
                    double r2 = sim.box.getLeafList().getAtom(i).getPosition().Mv1Squared(sim.box.getLeafList().getAtom(i<nPoints-1 ? i+1 : 0).getPosition());
                    x[0] = Math.sqrt(r2);
                    if (v!= 0) accHistTargRingy.putData(dataTarg);
                    dataRef.x = x[0];
                    accHistRefRingy.putData(dataRef);
                }
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        if (doHist) sim.integrator.getEventManager().addListener(histListenerRingy);

        if (false) {
            sim.box.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box); 
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            
            
            ColorScheme colorScheme = new ColorScheme() {
                
                public Color getAtomColor(IAtom a) {
                    float b=a.getLeafIndex()/((float)nPoints);
                    float r=1.0f-b;
                    return new Color(r, 0f, b);
                }
            };
            displayBox0.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.setAccumulatorBlockSize(1000);
            
            final JPanel panelParentGroup = new JPanel(new java.awt.GridBagLayout());
            GridBagConstraints gbc = new GridBagConstraints();
            gbc.gridx = 0;
            gbc.gridwidth = 2;
            final DisplayTextBox stepsBox = new DisplayTextBox();
            stepsBox.setLabel("steps");
            panelParentGroup.add(stepsBox.graphic(), gbc);

            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints);
            gbc.gridy = 1;
            panelParentGroup.add(jLabelPanelParentGroup, gbc);
            gbc.gridy = 2;
            gbc.gridwidth = 1;
            panelParentGroup.add(averageBox.graphic(), gbc);
            gbc.gridx = 1;
            panelParentGroup.add(errorBox.graphic(), gbc);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());
            
            IAction pushAnswer = new IAction() {
                public void actionPerformed() {
                    IData avgData = sim.accumulator.getData();
                    double ratio = ((DataGroup)avgData).getData(AccumulatorRatioAverageCovariance.RATIO.index).getValue(1);
                    double error = ((DataGroup)avgData).getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index).getValue(1);
                    data.x = ratio*refIntegral;
                    averageBox.putData(data);
                    data.x = error*Math.abs(refIntegral);
                    errorBox.putData(data);
                    
                    data.x = sim.integrator.getStepCount();
                    stepsBox.putData(data);
                }
                
                DataDouble data = new DataDouble();
            };
            IEtomicaDataInfo dataInfoSteps = new DataDouble.DataInfoDouble("B"+nPoints, Null.DIMENSION);
            stepsBox.putDataInfo(dataInfoSteps);
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, Null.DIMENSION);
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer, 1000));

            if (doHist) {
                DisplayPlot histPlot = new DisplayPlot();
                histPlot.setLabel("hist");
                simGraphic.add(histPlot);
                
                DataProcessor dpHistTarg = new DataProcessorR2();
                accHistTarg.addDataSink(dpHistTarg);
                dpHistTarg.setDataSink(histPlot.getDataSet().makeDataSink());
                histPlot.setLegend(new DataTag[]{dpHistTarg.getTag()}, "target");
                DataProcessor dpHistRef = new DataProcessorR2();
                accHistRef.addDataSink(dpHistRef);
                dpHistRef.setDataSink(histPlot.getDataSet().makeDataSink());
                histPlot.setLegend(new DataTag[]{dpHistRef.getTag()}, "ref");

                DataProcessor dpHistTargRingy = new DataProcessorR2();
                accHistTargRingy.addDataSink(dpHistTargRingy);
                dpHistTargRingy.setDataSink(histPlot.getDataSet().makeDataSink());
                histPlot.setLegend(new DataTag[]{dpHistTargRingy.getTag()}, "target (ring)");
                DataProcessor dpHistRefRingy = new DataProcessorR2();
                accHistRefRingy.addDataSink(dpHistRefRingy);
                dpHistRefRingy.setDataSink(histPlot.getDataSet().makeDataSink());
                histPlot.setLegend(new DataTag[]{dpHistRefRingy.getTag()}, "ref (ring)");

                histPlot.getPlot().setYLog(true);
            }

            return;
        }

        long t1 = System.currentTimeMillis();

        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        

        DataGroup allYourBase = (DataGroup)accumulator.getData();
        IData averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        IData errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        
        System.out.println();
        
        double avg = averageData.getValue(0);
        double err = errorData.getValue(0);
        System.out.print(String.format("target average: %20.15e error: %9.4e\n", avg, err));

        System.out.println();

        System.out.print(String.format("abs average: %20.15e  error: %9.4e\n", avg*refIntegral, err*Math.abs(refIntegral)));
            
//            double avg2 = averageData.getValue(1);
//            double err2 = errorData.getValue(1);
//            System.out.print(String.format("subset fraction of ref integral: %20.15e error: %9.4e\n", avg2, err2));

        
        if (doHist) {
            try {
                FileWriter fw = new FileWriter("hist"+nPoints+".dat");
                double[] x = histTarg.xValues();
                double[] y = histTarg.getHistogram();
                for (int i=0; i<x.length; i++) {
                    if (Double.isNaN(y[i]) || y[i] == 0) continue;
                    fw.write(x[i]+" "+y[i]+" "+y[i]/(x[i]*x[i])+"\n");
                }
                fw.close();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        System.out.println("time: "+(t2-t1)/1000.0);
        
        if (targetCluster instanceof ClusterWheatleyPartitionScreening) {
            ClusterWheatleyPartitionScreening tc = (ClusterWheatleyPartitionScreening)targetCluster;
            if(tc.sigCounter != null) {
                for (int i=3; i<6; i++) {
                    System.out.println("size "+i);
                    tc.sigCounter[i-1].print();
                }
                tc.sigCounter[nPoints-1].print();
                System.out.println("Fraction tabulated for each n");
                for(int i=0; i<nPoints; i++) {
                    System.out.println(tc.sigCounter[i].getn()+"\t"+tc.sigCounter[i].fractionTabulated()+"\t"+tc.sigCounter[i].getEntries());
                }
            }
        }
    }

    public static class DataProcessorR2 extends DataProcessor {
        DataFunction data;

        public DataPipe getDataCaster(IEtomicaDataInfo inDataInfo) {
            return null;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataInfoFunction("hist2", Null.DIMENSION, ((DataInfoFunction)inputDataInfo).getXDataSource());
            dataInfo.addTags(inputDataInfo.getTags());
            dataInfo.addTag(tag);
            data = new DataFunction(((DataInfoFunction)dataInfo).getArrayShape());
            return dataInfo;
        }

        protected IData processData(IData inputData) {
            double[] y = data.getData();
            IData xData = ((DataInfoFunction)dataInfo).getXDataSource().getIndependentData(0);
            for (int i=0; i<y.length; i++) {
                double r = xData.getValue(i);
                y[i] = inputData.getValue(i)/(r*r);
            }
            return data;
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 2;
        public long numSteps = 100000000;
        public static final int TREE = 0, CHAINS = 1, CHAIN_TREE = 5, RING_TREE = 7, RINGS = 8, RING_CHAIN_TREES = 9;
        public int ref = TREE;
        public boolean doHist = false;
        public double chainFrac = 1;
        public double ringFrac = 0.5;
        public boolean doResize = true;
        public String[] shapes = new String[0];
        public int[] nShapes = new int[0];
    }
    
}
