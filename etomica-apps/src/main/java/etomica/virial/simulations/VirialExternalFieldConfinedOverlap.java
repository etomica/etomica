/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.io.File;

import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graph.model.impl.MetadataImpl;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ClusterWeightSumWall;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirialExternalFieldOverlapConfined;
import etomica.virial.cluster.ExternalVirialDiagrams;
import etomica.virial.cluster.Standard;

/**
 * External Field simulation using Overlap2 sampling to evaluate cluster integrals
 */
public class VirialExternalFieldConfinedOverlap {


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
        double sigmaHSRef = 1;
        //
        double lambda = params.lambda;
        double epsilon = params.epsilon;
        //
        double [] walldistance = params.walldistance;
       

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
                
        System.out.println("External Field overlap sampling B(r1)"+nPoints+" at T="+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef); 
        //
        Potential2Spherical pTarget = new P2SquareWell(space, 1.0, lambda, epsilon, false);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        //then change meter and add input
        
        MetadataImpl.rootPointsSpecial=true;
        
        ExternalVirialDiagrams refDiagrams = new ExternalVirialDiagrams(nPoints, true, false, false);
        refDiagrams.setDoShortcut(true);
        ClusterSum refCluster = refDiagrams.makeRhoCluster(fRef, false);        
        
        refCluster.setTemperature(temperature);
        
        System.out.println(steps+" steps");
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);       
		
		//MeterVirialExternalFieldOverlapConfined meter = new MeterVirialExternalFieldOverlapConfined(refDiagrams, fRef, walldistance );
		MeterVirialExternalFieldOverlapConfined meter = new MeterVirialExternalFieldOverlapConfined(refDiagrams, fTarget, walldistance, temperature );
		
		final ClusterWeightSumWall targetSampleCluster =  new ClusterWeightSumWall(meter, nPoints);
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), temperature, new ClusterAbstract[] {refCluster,targetSampleCluster}, new ClusterWeight[] {sampleCluster, targetSampleCluster}, false);
        
        do {
      	  for (int i = 1; i<nPoints; i++){
      		  sim.box[1].getLeafList().get(i).getPosition().setX(2, sim.getRandom().nextDouble()-0.5);
      	  }
      	  sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
        } while (targetSampleCluster.value(sim.box[1]) == 0);
        
       
        DataProcessorFunction dividedbyPi = new DataProcessorFunction(null) {

			@Override
			protected IData processData(IData inputData) {
				double[] x = ((DataDoubleArray)inputData).getData();
				double[] y = ((DataDoubleArray)data).getData();
				double pi = targetSampleCluster.value(sim.box[1]);
				for(int i=0; i<x.length;i++){
					y[i] = x[i]/pi;
					
				}
				y[0] = sim.accumulators[1].getData(sim.accumulators[1].MOST_RECENT).getValue(1);
				return data;
			}
        	
		};
		
		DataPumpListener pump = new DataPumpListener(meter, dividedbyPi);
        sim.integrators[1].getEventManager().addListener(pump);
        AccumulatorAverageCovariance average = new AccumulatorAverageCovariance(true);
        dividedbyPi.setDataSink(average);      
             
        sim.initRefPref("RefPref", steps/40);        
        sim.equilibrate("RefPref", steps/20);
        
        System.out.println("equilibration finished");

        System.out.println("MC Move step sizes "+sim.mcMoveTranslate[0].getStepSize()+"  "+sim.mcMoveTranslate[1].getStepSize());
sim.getController().runActivityBlocking(new etomica.action.activity.ActivityIntegrate2(sim.integratorOS), steps);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());   
            
        sim.printResults(HSb[nPoints]);
        DataGroup allYourBase = (DataGroup)average.getData();
       
        for (int i=0; i < walldistance.length; i++){
            
        	IData covarianceData = allYourBase.getData(average.BLOCK_COVARIANCE.index);
            double correlationCoef = covarianceData.getValue(i+1)/Math.sqrt(covarianceData.getValue(0))/allYourBase.getData(average.STANDARD_DEVIATION.index).getValue(i+1);
            correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
            
            System.out.println(String.format("walldistance= %6.2f",walldistance[i]) + " target average: "+((DataDoubleArray)allYourBase.getData(average.AVERAGE.index)).getData()[i+1]
                                          
                                             +" error: "+((DataDoubleArray)allYourBase.getData(average.ERROR.index)).getData()[i+1]+" cov: "+correlationCoef);	
        }
       
    }

    

    /**
     * Inner class for parameters
     */
    public static class VirialExternalFieldParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 0.5;
        public long numSteps = 10000L;
        //
        public double epsilon = 1.0;
        public double lambda = 1.5;
        //
        public double [] walldistance = new double[]{5, 4, 3, 2, 1.1};
        
    }
}
