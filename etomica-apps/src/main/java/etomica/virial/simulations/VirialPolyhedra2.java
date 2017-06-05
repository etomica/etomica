package etomica.virial.simulations;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
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
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterChainHS;
import etomica.virial.ClusterSinglyConnected;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ClusterWeightUmbrella;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.MCMoveClusterPolyhedraChain;
import etomica.virial.MCMoveClusterPolyhedraTree;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralAtomic;
import etomica.virial.MeterVirialBD;
import etomica.virial.simulations.ShapeParser.ShapeData;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialPolyhedra2 {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.shape = "Cube";
            params.nPoints = 3;
            params.numSteps = 1000000L;
            params.ref = VirialHSParam.TREE;
            params.doHist = false;
            params.chainFrac = 0.5;
        }

        final String shape = params.shape;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;
        boolean doHist = params.doHist;
        final double chainFrac = params.chainFrac;

        System.out.println("Polyhedra ("+shape+") singly-connected sampling B"+nPoints);
        

        Space space = Space3D.getInstance();

        final P2SpheroPolyhedron p2 = new P2SpheroPolyhedron(space);

        SpeciesPolyhedron[] allSpecies = new SpeciesPolyhedron[1];
        ShapeData shapeData = ShapeParser.doParse("shape/"+shape+".dat", space);
        allSpecies[0] = new SpeciesPolyhedron(space, shapeData.vertices, 0.0, new ElementSimple("P"));
        double shsref = 2*((AtomTypeSpheroPolyhedron)allSpecies[0].getAtomType(0)).getOuterRadius();
        final double sigmaHSRef = shsref;
        double B2 = shapeData.B2;

        if (nPoints==2) {
            System.out.println("abs value: "+B2);
            System.exit(0);
        }

        final double[][] uValues = new double[nPoints][nPoints];
        IPotentialAtomic p2Wrapper = new IPotentialAtomic() {
            
            public void setBox(Box box) {
                p2.setBox(box);
            }

            public int nBody() {
                return 2;
            }
            
            public double getRange() {
                return sigmaHSRef;
            }
            
            public double energy(IAtomList atoms) {
                int i0 = atoms.getAtom(0).getLeafIndex();
                int i1 = atoms.getAtom(1).getLeafIndex();
                if (Double.isNaN(uValues[i0][i1])) {
                    double u = p2.energy(atoms);
                    uValues[i0][i1] = uValues[i1][i0] = u;
                }
                return uValues[i0][i1];
            }
        };
        MayerFunction fTarget = new MayerGeneralAtomic(p2Wrapper);

        ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fTarget);
        
        System.out.println("reference diameter: "+sigmaHSRef);

        targetCluster.setTemperature(1.0);

        ClusterAbstract refCluster = null;
        long numDiagrams = 0;
        
        double ri = 0;
        if (ref == VirialHSParam.TREE) {
            System.out.println("using a tree reference");
            refCluster = new ClusterSinglyConnected(nPoints, fTarget);
            numDiagrams = ((ClusterSinglyConnected)refCluster).numDiagrams();
            ri = numDiagrams*Math.pow(2*B2, nPoints-1);
        }
        else if (ref == VirialHSParam.CHAINS) {
            System.out.println("using a chain reference");
            refCluster = new ClusterChainHS(nPoints, fTarget);
            numDiagrams = ((ClusterChainHS)refCluster).numDiagrams();
            ri = numDiagrams*Math.pow(2*B2, nPoints-1);
        }
        else if (ref == VirialHSParam.CHAIN_TREE) {
            System.out.println("using a chain/tree reference ("+chainFrac+" chains)");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fTarget);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fTarget);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }
            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{chainFrac/(SpecialFunctions.factorial(nPoints)/2),(1-chainFrac)/numTreeDiagrams});
            ri = Math.pow(-2*B2, nPoints-1);
        }


        // (nPoints-1)! is simply not included by ClusterWheatley, so do that here.
        final double refIntegral = ri;
        System.out.println("reference integral: "+refIntegral);
        refCluster.setTemperature(1.0);



        System.out.println(steps+" steps");

        final SimulationVirial sim = new SimulationVirial(space, allSpecies, new int[]{nPoints}, 1.0,ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, new ClusterAbstract[]{targetCluster});
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
            MCMoveClusterPolyhedraTree mcMoveTree = new MCMoveClusterPolyhedraTree(sim.getRandom(), space, sigmaHSRef, p2, uValues);
            sim.integrator.getMoveManager().addMCMove(mcMoveTree);
        }
        else if (ref == VirialHSParam.CHAINS) {
            MCMoveClusterPolyhedraChain mcMoveChain = new MCMoveClusterPolyhedraChain(sim.getRandom(), space, sigmaHSRef, p2, uValues);
            sim.integrator.getMoveManager().addMCMove(mcMoveChain);
        }
        else if (ref == VirialHSParam.CHAIN_TREE) {
            MCMoveClusterPolyhedraTree mcMoveTree = new MCMoveClusterPolyhedraTree(sim.getRandom(), space, sigmaHSRef, p2, uValues);
            sim.integrator.getMoveManager().addMCMove(mcMoveTree);
            sim.integrator.getMoveManager().setFrequency(mcMoveTree, 1-chainFrac);
            MCMoveClusterPolyhedraChain mcMoveChain = new MCMoveClusterPolyhedraChain(sim.getRandom(), space, sigmaHSRef, p2, uValues);
            sim.integrator.getMoveManager().addMCMove(mcMoveChain);
            sim.integrator.getMoveManager().setFrequency(mcMoveChain, chainFrac);
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
        IIntegratorListener histListener = new IIntegratorListener() {
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
        IIntegratorListener histListenerRingy = new IIntegratorListener() {
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
                    double ratio = ((DataGroup)avgData).getData(sim.accumulator.RATIO.index).getValue(1);
                    double error = ((DataGroup)avgData).getData(sim.accumulator.RATIO_ERROR.index).getValue(1);
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
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);
        
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
        public static final int TREE = 0, CHAINS = 1, CHAIN_TREE = 2;
        public int ref = TREE;
        public boolean doHist = false;
        public double chainFrac = 1;
        public String shape = "";
    }
    
}
