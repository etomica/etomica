package etomica.virial.simulations;


import java.awt.Color;

import etomica.action.IAction;
import etomica.api.IPotentialMolecular;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterCoupledFlipped;
import etomica.virial.ClusterSumPolarizable;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryWaterGCPM;
import etomica.virial.cluster.Standard;


public class VirialGCPM {

    public static void main(String[] args) {

        VirialGCPMParam params = new VirialGCPMParam();
        if (args.length > 0) {
            ReadParameters paramReader = new ReadParameters(args[0], params);
            paramReader.readParameters();
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        int numSubSteps = 1000;
        double deltaCut = 100;

        double sigmaHSRef = 3.2;
        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Water overlap sampling B"+nPoints+" at T="+temperature);
        temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println(steps+" steps ("+steps/1000+" blocks of 1000)");
        steps /= 1000;

        Space space = Space3D.getInstance();

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        final IPotentialMolecular pTarget = new PNWaterGCPM(space);
        
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        ClusterAbstract targetCluster = Standard.virialClusterPolarizable(nPoints, fTarget, nPoints>3, false);
        ((ClusterSumPolarizable)targetCluster).setDeltaCut(deltaCut);
        targetCluster = new ClusterCoupledFlipped(targetCluster, space);

   	    ClusterWeight sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster);

        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster);

       
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refSample.setTemperature(temperature);


        final SimulationVirialOverlapRejected sim = new SimulationVirialOverlapRejected(space,new SpeciesFactoryWaterGCPM(), 
                temperature, new ClusterAbstract[]{refCluster,targetCluster},new ClusterWeight[]{refSample,sampleCluster1}, false);
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesWater4P species = (SpeciesWater4P)sim.getSpeciesManager().getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null,20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize());
        
//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

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
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
	}



    /**
     * Inner class for parameters
     */
    public static class VirialGCPMParam extends ParameterBase {
        public int nPoints = 5;
        public double temperature = 350;   // Kelvin
        public long numSteps = 100000;
    }
}
