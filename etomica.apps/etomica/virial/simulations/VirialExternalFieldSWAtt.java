package etomica.virial.simulations;

import java.io.File;

import etomica.api.IAtomType;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graph.model.impl.MetadataImpl;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumExternalField;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirialExternalField;
import etomica.virial.MeterVirialExternalFieldConfined;
import etomica.virial.MeterVirialExternalFieldSW;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.ExternalVirialDiagrams;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.VirialExternalField3.VirialExternalFieldParam;

/**
 * External Field simulation using Direct sampling to evaluate cluster integrals
 */
public class VirialExternalFieldSWAtt {


	public static void main(String[] args) {
    	MetadataImpl.rootPointsSpecial=true;
    	VirialExternalFieldParam params = new VirialExternalFieldParam();
    	if (args.length > 0) {
	        if (new File(args[0]).exists()) {
	           
	            ReadParameters readParameters = new ReadParameters(args[0], params);
	            readParameters.readParameters();
	            args = (String[])Arrays.removeObject(args, args[0]);
	        }
	        if (args.length > 0) {
	        	ParseArgs parseArgs = new ParseArgs(params);
	        	parseArgs.parseArgs(args);
	        }
    	}        
        runVirial(params);
    }
    
    public static void runVirial(VirialExternalFieldParam params) {
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double lambda = params.lambda;
        double lambdaWF = params.lambdaWF;
        double epsilon = params.epsilon;
        double epsilonWF = params.epsilonWF;
        
        
        final double[] SWb = new double[6];
        SWb[2] = -1.0*4*Math.PI/6*Math.pow(1, 3)*(1-(Math.pow(lambda, 3)-1)*(Math.exp(epsilon/temperature)-1));
        SWb[3] = 2.0*Math.pow(SWb[2], 2)-0.5*2*Math.pow(Math.PI/6, 2)*(5-(Math.pow(lambda, 6)-18*Math.pow(lambda, 4)+32*Math.pow(lambda, 3)-15)*(Math.exp(epsilon/temperature)-1)+(-2*Math.pow(lambda, 6)+36*Math.pow(lambda, 4)-32*Math.pow(lambda, 3)-18*Math.pow(lambda, 2)+16)*Math.pow(Math.exp(epsilon/temperature)-1, 2)-(6*Math.pow(lambda, 6)-18*Math.pow(lambda, 4)+18*Math.pow(lambda, 2)-6)*Math.pow(Math.exp(epsilon/temperature)-1, 3));
        //SWb[4] = 2*Math.pow(Math.PI/6, 2)*(5-(Math.pow(lambda, 6)-18*Math.pow(lambda, 4)+32*Math.pow(lambda, 3)-15)*(Math.exp(1/temperature)-1)+(-2*Math.pow(lambda, 6)+36*Math.pow(lambda, 4)-32*Math.pow(lambda, 3)-18*Math.pow(lambda, 2)+16)*Math.pow(Math.exp(1/temperature)-1, 2)-(6*Math.pow(lambda, 6)-18*Math.pow(lambda, 4)+18*Math.pow(lambda, 2)-6)*Math.pow(Math.exp(1/temperature)-1, 3));
        //SWb[3] = 2.0*Math.pow(SWb[2], 2)-0.5*Math.pow(2*Math.PI/3*6.0221415*Math.pow(10, 23), 2)/8*(5-17*(Math.exp(1/temperature)-1)-Math.pow(Math.exp(1/temperature)-1, 2)*(48+18*Math.pow(lambda, 2)-32*Math.pow(lambda, 3))-Math.pow(Math.exp(1/temperature)-1, 3)*(5*Math.pow(lambda, 6)-32*Math.pow(lambda, 3)+18*Math.pow(lambda, 2)+26));
       
        System.out.println("b2SW: "+SWb[2]);
        System.out.println("b3SW: "+SWb[3]); 
        //System.out.println("B3SW: "+SWb[4]); 
                        
        System.out.println("External Field direct sampling b"+nPoints+" at T="+temperature);
        System.out.println("lambda: "+lambda);
        System.out.println("lambdaWF: "+lambdaWF);
		
        Space space = Space3D.getInstance();
        
        Potential2Spherical pSW = new P2SquareWell(space,1.0,lambda,epsilon, false);
        MayerGeneralSpherical fRef = new MayerGeneralSpherical(pSW);
                             
        ExternalVirialDiagrams refDiagrams = new ExternalVirialDiagrams(nPoints, false, false);
        refDiagrams.setDoShortcut(true);
        ClusterSum refCluster = refDiagrams.makeRhoCluster(fRef, false);
        
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps");
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        final SimulationVirial sim = new SimulationVirial(space,new SpeciesFactorySpheres(), temperature, sampleCluster, refCluster,new ClusterAbstract[0]);
        double[] wallposition = new double[(int)Math.round((lambdaWF+(nPoints-1)*lambda-0.5)*100)+1];
        for (int i=0; i < wallposition.length; i++){
        	wallposition[i] = -lambdaWF-lambda*(nPoints-1)+0.01*i;
        }
        //double temperatureWF = temperature/epsilonWF;
        MeterVirialExternalFieldSW meter = new MeterVirialExternalFieldSW(refCluster, wallposition, lambdaWF, temperature, epsilonWF);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        sim.setAccumulator(new AccumulatorRatioAverageCovariance(steps/100));
        
      
       
      
        sim.equilibrate(steps/40);
        
        System.out.println("equilibration finished");



        sim.ai.setMaxSteps(steps);
        
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
        
        sim.getController().actionPerformed();

        
       
      
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        double sum=0;
        for (int i=0; i < wallposition.length; i++){
            sum+=(((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[i]*SWb[nPoints]-SWb[nPoints])*0.01;
            System.out.println(String.format("wallposition= %6.2f",wallposition[i]) + " ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[i]
                                             +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[i] + " reference average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[i]
                                             +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[i]
                                             +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[i]);	
        }
        System.out.println("sum="+sum);
        System.out.println("surfacevirial "+" ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[wallposition.length]
                           +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[wallposition.length] + " reference average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[wallposition.length]
                           +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[wallposition.length]
                           +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[wallposition.length]);
        System.out.println("b"+nPoints+"="+((((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[wallposition.length])*SWb[nPoints]-(lambdaWF+lambda*(nPoints-1)-0.5)*SWb[nPoints])
                         +" error: "+SWb[nPoints]*((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[wallposition.length]);
    }
	

    /**
     * Inner class for parameters
     */
    public static class VirialExternalFieldParam extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 30;
        public long numSteps = 10000000L;
        public double lambda = 1.7;
        public double lambdaWF = 1.1;
        public double epsilon = 1.0;
        public double epsilonWF = 2.0;
        
       
    }
}