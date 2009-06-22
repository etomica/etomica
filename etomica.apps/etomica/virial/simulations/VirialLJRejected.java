package etomica.virial.simulations;

import java.awt.Color;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IData;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Null;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.Volume;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.util.Constants.CompassDirection;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJRejected {


    public static void main(String[] args) {
        VirialLJParam params = new VirialLJParam();
        if (args.length > 0) {
            params.writeRefPref = true;
            ReadParameters readParameters = new ReadParameters(args[0], params);
            readParameters.readParameters();
        }
        
        runVirial(params);
    }
    
    public static void runVirial(VirialLJParam params) {
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
        System.out.println("B8HS: "+HSB[8]+" = 0.004164 B2HS^7");
        System.out.println("Lennard Jones overlap sampling B"+nPoints+" at T="+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        MayerESpherical eTarget = new MayerESpherical(pTarget);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		
        final SimulationVirialOverlapRejected sim = new SimulationVirialOverlapRejected(space,new SpeciesFactorySpheres(), temperature,refCluster,targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(steps*10);
        
        if (true) {
            double size = 5;
            sim.box[0].getBoundary().setDimensions(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setDimensions(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box[0]).canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box[1]).canvas).setBackgroundColor(Color.WHITE);
            
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            simGraphic.getDisplayBox(sim.box[0]).setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            simGraphic.getDisplayBox(sim.box[1]).setColorScheme(colorScheme);
            
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 1000);
                    sim.equilibrate(null, 2000);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            final DisplayTextBox averageBox = new DisplayTextBox("Average", Quantity.SIM_UNIT);
            final DisplayTextBox errorBox = new DisplayTextBox("Error", Quantity.SIM_UNIT);
            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
            final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());
            
            IAction pushAnswer = new IAction() {
                public void actionPerformed() {
                    double ratio = sim.dsvo.getDataAsScalar() * HSB[nPoints];
                    double error = sim.dsvo.getError() * HSB[nPoints];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
                
                DataDouble data = new DataDouble();
            };
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
            Unit unit = Null.UNIT;
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            averageBox.setUnit(unit);
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            errorBox.setUnit(unit);
            IntegratorListenerAction pushAnswerListener = new IntegratorListenerAction(pushAnswer);
            pushAnswerListener.setInterval(1);
            sim.integratorOS.getEventManager().addListener(pushAnswerListener);
            
            return;
        }

        
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
        
        System.out.println("equilibration finished");

        IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double ratio = sim.dsvo.getDataAsScalar();
                double error = sim.dsvo.getError();
                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval((int)(steps/10));
        sim.integratorOS.getEventManager().addListener(progressReportListener);
        
//        VirialHistogram virialHistogram = new VirialHistogram(new P2HardSphere(space, sigmaHSRef, false), pTarget, sim.box[1]);
//        virialHistogram.setBinFac(100);
//        sim.integrators[1].addIntervalAction(virialHistogram);

        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().actionPerformed();
        
//        long[][] histogram = virialHistogram.getHistogram();
//        int numNegBins = virialHistogram.getNumNegBins();
//        double zeroOffset = virialHistogram.getZeroOffset();
//        double binFac = virialHistogram.getBinFac();
//        if (false) {
//        try {
//            FileWriter writer = new FileWriter("nover.dat");
//            for (int i=0; i<numNegBins; i++) {
//                if (histogram[0][i] == 0) continue;
//                writer.write(-(Math.exp(i/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[0][i]+"\n");
//            }
//            for (int i=numNegBins; i<histogram[0].length; i++) {
//                if (histogram[0][i] == 0) continue;
//                writer.write((Math.exp((i-numNegBins)/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[0][i]+"\n");
//            }
//            writer.close();
//            writer = new FileWriter("over.dat");
//            for (int i=0; i<numNegBins; i++) {
//                if (histogram[1][i] == 0) continue;
//                writer.write(-(Math.exp(i/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[1][i]+"\n");
//            }
//            for (int i=numNegBins; i<histogram[1].length; i++) {
//                if (histogram[1][i] == 0) continue;
//                writer.write((Math.exp((i-numNegBins)/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[1][i]+"\n");
//            }
//            writer.close();
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
//        }

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());
        
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
        IData ratioData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index);
        IData ratioErrorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index);
        IData averageData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        IData stdevData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.ERROR.index);
        System.out.println("hard sphere ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("hard sphere   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("hard sphere overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.ERROR.index);
        System.out.println("lennard jones ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("lennard jones average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("lennard jones overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 1.3;
        public long numSteps = 10000;
        public double sigmaHSRef = 1.5;
        public boolean writeRefPref = false;
    }
}
