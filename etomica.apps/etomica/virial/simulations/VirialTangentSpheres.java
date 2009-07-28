package etomica.virial.simulations;

import etomica.api.ISpecies;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntergroup;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.potential.P2Exp6Buckingham;
import etomica.potential.P2HardSphere;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.cluster.Standard;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialTangentSpheres {


    public static void main(String[] args) {
        VirialTangentSpheresParam params = new VirialTangentSpheresParam();
        if (args.length > 0) {
            ReadParameters paramReader = new ReadParameters(args[0], params);
            paramReader.readParameters();
        }
        int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        String model = params.model;
        double bondL = params.bondL;
        double sigmaHSRef = Math.pow(nSpheres,1.0/3.0);
        if (temperature > 0) {
            sigmaHSRef += 0.5;
        }
        double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        PotentialGroup pTargetGroup = new PotentialGroup(2);
        Potential2 p2 = null;
        if (model.equals("hard")) {
            if (temperature == 0) {
                System.out.println(nSpheres+"-mer hard chain B"+nPoints);
                p2 = new P2HardSphere(space, 1.0, false);
                temperature = 1.0;
            }
            else {
                System.out.println(nSpheres+"-mer sqw chains B"+nPoints+" at "+temperature);
                p2 = new P2SquareWell(space, 1.0, 1.5, 1.0, false);
            }
        }
        else if (model.equals("LJ")) {
            System.out.println(nSpheres+"-mer LJ chains B"+nPoints+" at "+temperature);
            p2 = new P2LennardJones(space, 1.0, 1.0);
        }
        else if (model.equals("exp6")) {
            double alpha = params.exp6Alpha;
            System.out.println(nSpheres+"-mer exp6(alpha="+alpha+") chains B"+nPoints+" at "+temperature);
            double lambda = 0;
            double rm = 0;
            if (alpha == 15) {
                lambda = 0.1682455;
                rm = 1.0/0.894170;
            }
            else if (alpha == 30) {
                lambda = 0.1465604;
                rm = 1.0 / 0.932341;
            }
            else if (alpha == 100) {
                lambda = 6.248806e-7;
                rm = 1.0/0.970041;
            }
            p2 = new P2Exp6Buckingham(space, 1.0, alpha, rm, lambda*rm);
        }
        else {
            throw new RuntimeException("Unknown model "+model);
        }
        System.out.println("bond length: "+bondL);
        pTargetGroup.addPotential(p2, new ApiIntergroup());
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        SimulationVirialOverlapRejected sim = new SimulationVirialOverlapRejected(space,new SpeciesFactoryTangentSpheres(nSpheres,
                new ConformationLinear(space,bondL)), temperature,refCluster,targetCluster, true);
        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(10*steps);
        
        if (nSpheres > 2) {
            PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
            pIntra.addPotential(p2,ApiBuilder.makeNonAdjacentPairIterator());
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.species});
        }
        
        if (p2 instanceof P2HardSphere) {
            temperature = 0;
        }
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle[0]==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle[1]==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        
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
        System.out.println("reference average: "+averageData.getValue(0)
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
    public static class VirialTangentSpheresParam extends ParameterBase {
        public int nPoints = 2;
        public int nSpheres = 2;
        public double temperature = 500.0/114.0;
        public long numSteps = 10000;
        public String model = "LJ";
        public double exp6Alpha = 15;
        public double bondL = 1.54/3.93;
    }
}
