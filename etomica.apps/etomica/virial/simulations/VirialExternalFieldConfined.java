package etomica.virial.simulations;

import java.io.File;

import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graph.model.impl.MetadataImpl;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirialExternalFieldConfined;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.ExternalVirialDiagrams;
import etomica.virial.cluster.Standard;

/**
 * External Field simulation using Direct sampling to evaluate cluster integrals
 */
public class VirialExternalFieldConfined {


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
        double [] walldistance = params.walldistance;
        long steps = params.numSteps;
        double sigmaHSRef = 1;
       

        final double[] HSb = new double[9];
        HSb[2] = -1.0*Standard.B2HS(sigmaHSRef);
        HSb[3] = 2.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)-0.5*Standard.B3HS(sigmaHSRef);
        HSb[4] = 1.0/3.0*(-16.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)+9.0*Standard.B2HS(sigmaHSRef)*Standard.B3HS(sigmaHSRef)-Standard.B4HS(sigmaHSRef));
        HSb[5] = 1.0/24.0*(400.0*Math.pow(Standard.B2HS(sigmaHSRef), 4)-360.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B3HS(sigmaHSRef)+27.0*Math.pow(Standard.B3HS(sigmaHSRef), 2)+64.0*Standard.B2HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)-6*Standard.B5HS(sigmaHSRef));       
        HSb[6] = 1.0/10.0*(-576.0*Math.pow(Standard.B2HS(sigmaHSRef), 5)+720.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)*Standard.B3HS(sigmaHSRef)-135.0*Standard.B2HS(sigmaHSRef)*Math.pow(Standard.B3HS(sigmaHSRef), 2)-160.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B4HS(sigmaHSRef)+20.0*Standard.B3HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)+25.0*Standard.B2HS(sigmaHSRef)*Standard.B5HS(sigmaHSRef)-2.0*Standard.B6HS(sigmaHSRef));
        HSb[7] = 1.0/720.0*(153664.0*Math.pow(Standard.B2HS(sigmaHSRef), 6)-246960.0*Math.pow(Standard.B2HS(sigmaHSRef), 4)*Standard.B3HS(sigmaHSRef)+79380.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Math.pow(Standard.B3HS(sigmaHSRef), 2)-2835.0*Math.pow(Standard.B3HS(sigmaHSRef), 3)+62720.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)*Standard.B4HS(sigmaHSRef)-20160.0*Standard.B2HS(sigmaHSRef)*Standard.B3HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)+640.0*Math.pow(Standard.B4HS(sigmaHSRef), 2)-12600.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B5HS(sigmaHSRef)+1350.0*Standard.B3HS(sigmaHSRef)*Standard.B5HS(sigmaHSRef)+1728.0*Standard.B2HS(sigmaHSRef)*Standard.B6HS(sigmaHSRef)-120.0*Standard.B7HS(sigmaHSRef));
        HSb[8] = 1.0/315.0*(-262144.0*Math.pow(Standard.B2HS(sigmaHSRef), 7)+516096.0*Math.pow(Standard.B2HS(sigmaHSRef), 5)*Standard.B3HS(sigmaHSRef)-241920.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)*Math.pow(Standard.B3HS(sigmaHSRef), 2)+22680.0*Standard.B2HS(sigmaHSRef)*Math.pow(Standard.B3HS(sigmaHSRef), 3)-143360.0*Math.pow(Standard.B2HS(sigmaHSRef), 4)*Standard.B4HS(sigmaHSRef)+80640.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B3HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)-3780.0*Math.pow(Standard.B3HS(sigmaHSRef), 2)*Standard.B4HS(sigmaHSRef)-4480.0*Math.pow(Standard.B4HS(sigmaHSRef), 2)*Standard.B2HS(sigmaHSRef)+33600.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)*Standard.B5HS(sigmaHSRef)-9450.0*Standard.B2HS(sigmaHSRef)*Standard.B3HS(sigmaHSRef)*Standard.B5HS(sigmaHSRef)+525.0*Standard.B4HS(sigmaHSRef)*Standard.B5HS(sigmaHSRef)-6048.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B6HS(sigmaHSRef)+567.0*Standard.B3HS(sigmaHSRef)*Standard.B6HS(sigmaHSRef)+735.0*Standard.B2HS(sigmaHSRef)*Standard.B7HS(sigmaHSRef)-45.0*Standard.B8HS(sigmaHSRef));
        
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("b2HS: "+HSb[2]);
        System.out.println("b3HS: "+HSb[3]);
        System.out.println("b4HS: "+HSb[4]);
        System.out.println("b5HS: "+HSb[5]);
        System.out.println("b6HS: "+HSb[6]);
        System.out.println("b7HS: "+HSb[7]);
        System.out.println("b8HS: "+HSb[8]);
                
        System.out.println("External Field direct sampling b"+nPoints+" at T="+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);                               
       
     
        ExternalVirialDiagrams refDiagrams = new ExternalVirialDiagrams(nPoints, false, false);
        refDiagrams.setDoShortcut(true);
        ClusterSum refCluster = refDiagrams.makeRhoCluster(fRef, false);
        
        
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps");
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        final SimulationVirial sim = new SimulationVirial(space,new SpeciesFactorySpheres(), temperature,sampleCluster, refCluster,new ClusterAbstract[0]);
        double[] wallposition = new double[(nPoints-1)*100+1];
        for (int i=0; i < wallposition.length; i++){
        	wallposition[i] = 0.5-nPoints+0.01*i;
        }
        MeterVirialExternalFieldConfined meter = new MeterVirialExternalFieldConfined(refCluster, wallposition, walldistance);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        sim.setAccumulator(new AccumulatorRatioAverageCovariance(steps/100));
        
      
       
      
        sim.equilibrate(steps/40);
        
        System.out.println("equilibration finished");



        sim.ai.setMaxSteps(steps);
        
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
        
        sim.getController().actionPerformed();

        
       
      
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        int n=walldistance.length+1;
        for (int j=0; j<walldistance.length+1; j++){
	        double sum=0;
	        for (int i=0; i < wallposition.length; i++){
	            sum+=(((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[i*n+j+1]*HSb[nPoints]-HSb[nPoints])*0.01;
	            System.out.println(String.format("wallposition= %6.2f",wallposition[i]) + " ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[i*n+j+1]
	                                             +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index)).getData()[i*n+j+1] + " reference average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.AVERAGE.index)).getData()[i*n+j+1]
	                                             +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.STANDARD_DEVIATION.index)).getData()[i*n+j+1]
	                                             +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.ERROR.index)).getData()[i*n+j+1]);	
	        }
	      
	        System.out.println(String.format("walldistance= %6.2f", walldistance[j]) + " surfacevirial "+" ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[n*wallposition.length+j]
	                           +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index)).getData()[n*wallposition.length+j] + " reference average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.AVERAGE.index)).getData()[n*wallposition.length+j]
	                           +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.STANDARD_DEVIATION.index)).getData()[n*wallposition.length+j]
	                           +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.ERROR.index)).getData()[n*wallposition.length+j]);
	       
	    }
    }	

    /**
     * Inner class for parameters
     */
    public static class VirialExternalFieldParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 1.5;
        public double [] walldistance = new double[]{3.5,6,8,10};
        public long numSteps = 1000000L;
        
    }
}
