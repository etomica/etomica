package etomica.virial.simulations;

import java.awt.Color;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.model.Metadata;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.AllIsomorphs;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.AllIsomorphs.AllIsomorphsParameters;
import etomica.graph.property.IsFFT;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.Constants.CompassDirection;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterChainWheatley;
import etomica.virial.ClusterDifference;
import etomica.virial.ClusterPY;
import etomica.virial.ClusterSinglyConnected;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ClusterWheatleyPartitionScreening;
import etomica.virial.MCMoveClusterAtomHSChain;
import etomica.virial.MCMoveClusterAtomHSTree;
import etomica.virial.MayerHardSphere;
import etomica.virial.PYGenerator;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialHS {


    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        ParseArgs.doParseArgs(params, args);
        
        runVirial(params);
    }
    
    public static void runVirial(VirialHSParam params) {
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        int ref = params.ref;
        int subtractGraphs = params.subtractGraphs;
        double sigmaHS = 1.0;

        double litHSB = Double.NaN;
        try {
            litHSB = Standard.BHS(nPoints, sigmaHS);
        }
        catch (RuntimeException e) {}
        double bhs = -(4.0/3.0)*Math.PI*sigmaHS*sigmaHS*sigmaHS;
        final double refIntegral = Math.pow(bhs, nPoints-1);

        System.out.println("HS singly-connected sampling B"+nPoints);
		
        Space space = Space3D.getInstance();
        
        VirialDiagrams diagrams = new VirialDiagrams(nPoints, false, false);
        diagrams.setAllPermutations(true);
        diagrams.setDoReeHoover(false);
        diagrams.setDoShortcut(true);
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHS);
        ClusterAbstract targetCluster = new ClusterWheatleyPartitionScreening(nPoints, fRef);
        Set<Graph> myGraphs = new HashSet<Graph>();
        if (subtractGraphs == VirialHSParam.PY || subtractGraphs == VirialHSParam.ICPY) {
            System.out.println("Computing "+(subtractGraphs==VirialHSParam.ICPY ? "incremental":"")+" correction to PY");

            MetadataImpl.rootPointsSpecial = true;
            Set<Graph> icGraphs = (subtractGraphs == VirialHSParam.ICPY) ? PYGenerator.getICPYCorrection((byte)nPoints, true) : PYGenerator.getPYCorrection((byte)nPoints);
            MetadataImpl.rootPointsSpecial = false;
            for (Graph g: icGraphs) {
                g.getNode((byte)0).setType(Metadata.TYPE_NODE_FIELD);
                g.getNode((byte)1).setType(Metadata.TYPE_NODE_FIELD);
            }
            IsoFree isoFree = new IsoFree();
            AllIsomorphs allIso = new AllIsomorphs();
            icGraphs = isoFree.apply(icGraphs, null);
            icGraphs = allIso.apply(icGraphs, new AllIsomorphsParameters(true));
            MulScalarParameters msp = new MulScalarParameters(1, -nPoints);
            MulScalar mulScalar = new MulScalar();
            myGraphs = mulScalar.apply(icGraphs, msp);
            diagrams.makeVirialDiagrams();
            targetCluster = diagrams.makeVirialCluster(myGraphs, fRef, null);
        }
        else if (subtractGraphs == VirialHSParam.PY2) {
            System.out.println("Computing PYC recusively");
            targetCluster = new ClusterPY(nPoints, fRef);
        }
        else if (subtractGraphs == VirialHSParam.FFT) {
            System.out.println("Subtracting FFT diagrams");
            if (subtractGraphs == 1) throw new RuntimeException("nope");
            IsFFT isFFT = new IsFFT();
            diagrams.makeVirialDiagrams();
            for (Graph g : diagrams.getMSMCGraphs(true, false)) {
                if (subtractGraphs == 1 && g.edgeCount() == nPoints) {
                    myGraphs.add(g);
                }
                else if (subtractGraphs == 2 && isFFT.check(g)) {
                    myGraphs.add(g);
                }
            }
            targetCluster = new ClusterDifference(targetCluster, new ClusterAbstract[]{diagrams.makeVirialCluster(myGraphs, fRef, null)});
        }
        targetCluster.setTemperature(1.0);
        
        ClusterAbstract refCluster = null;
        long numDiagrams = 0;
        if (ref == VirialHSParam.TREE) {
            System.out.println("using a tree reference");
            refCluster = new ClusterSinglyConnected(nPoints, fRef);
            numDiagrams = ((ClusterSinglyConnected)refCluster).numDiagrams();
        }
        else {
            System.out.println("using a chain reference");
            refCluster = new ClusterChainWheatley(nPoints, fRef);
            numDiagrams = ((ClusterChainWheatley)refCluster).numDiagrams();
        }
        refCluster.setTemperature(1.0);

        System.out.println(steps+" steps");
		
        final SimulationVirial sim = new SimulationVirial(space,new SpeciesSpheresMono(space, new ElementSimple("A")), 1.0,ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, new ClusterAbstract[]{targetCluster});

        if (false) {
            sim.box.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box); 
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            
            
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.setAccumulatorBlockSize(1000);
            
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints);
            final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
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
                }
                
                DataDouble data = new DataDouble();
            };
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));

            return;
        }
        
        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        if (ref == VirialHSParam.TREE) {
            MCMoveClusterAtomHSTree mcMoveHS = new MCMoveClusterAtomHSTree(sim.getRandom(), space, 1);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        }
        else {
            MCMoveClusterAtomHSChain mcMoveHS = new MCMoveClusterAtomHSChain(sim.getRandom(), space, 1);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        }

        long t1 = System.currentTimeMillis();
        
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        
        if (!Double.isNaN(litHSB) && subtractGraphs == VirialHSParam.NONE) System.out.println("lit value "+litHSB);
        
        sim.printResults(refIntegral*numDiagrams);
        System.out.println("time: "+(t2-t1)/1000.0);
        
        ClusterWheatleyPartitionScreening tc = (ClusterWheatleyPartitionScreening)targetCluster;
        if(tc.sigCounter != null) {
            tc.sigCounter[nPoints-1].print();
            System.out.println("Fraction tabulated for each n");
            for(int i=0; i<nPoints; i++) {
                System.out.println(tc.sigCounter[i].getn()+"\t"+tc.sigCounter[i].fractionTabulated()+"\t"+tc.sigCounter[i].getEntries());
            }
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 8;
        public long numSteps = 1000000;
        public int ref = 1;
        public static final int TREE = 0, CHAINS = 1;
        public int subtractGraphs = 0;
        public static final int NONE = 0, RINGS = 1, FFT = 2, PY = 3, ICPY = 4, PY2 = 5;
    }
}
