package etomica.virial.simulations;

import java.io.File;

import etomica.api.IAtomType;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graph.model.impl.MetadataImpl;
import etomica.potential.P1HardBoundary;
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
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.ExternalVirialDiagrams;
import etomica.virial.cluster.Standard;

/**
 * External Field simulation using Direct sampling to evaluate cluster integrals
 */
public class VirialExternalField {


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
        double z0=params.z0;

        final double[] HSb = new double[9];
        HSb[2] = -1.0*Standard.B2HS(sigmaHSRef);
        HSb[3] = 2.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)-0.5*Standard.B3HS(sigmaHSRef);
        HSb[4] = 1.0/3.0*(-16.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)+9.0*Standard.B2HS(sigmaHSRef)*Standard.B3HS(sigmaHSRef)-Standard.B4HS(sigmaHSRef));
        HSb[5] = 1.0/24.0*(400.0*Math.pow(Standard.B2HS(sigmaHSRef), 4)-360.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B3HS(sigmaHSRef)+27.0*Math.pow(Standard.B3HS(sigmaHSRef), 2)+64.0*Standard.B2HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)-6*Standard.B5HS(sigmaHSRef));       
        HSb[6] = 1.0/10.0*(-576.0*Math.pow(Standard.B2HS(sigmaHSRef), 5)+720.0*Math.pow(Standard.B2HS(sigmaHSRef), 3)*Standard.B3HS(sigmaHSRef)-135.0*Standard.B2HS(sigmaHSRef)*Math.pow(Standard.B3HS(sigmaHSRef), 2)-160.0*Math.pow(Standard.B2HS(sigmaHSRef), 2)*Standard.B4HS(sigmaHSRef)+20.0*Standard.B3HS(sigmaHSRef)*Standard.B4HS(sigmaHSRef)+25.0*Standard.B2HS(sigmaHSRef)*Standard.B6HS(sigmaHSRef)-2.0*Standard.B6HS(sigmaHSRef));
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
       
        /*int [][][]bondlist=new int [][][]{{{0,1}}};
        ClusterBonds clusterbonds=new ClusterBonds(2,bondlist);
        ClusterSumExternalField targetCluster = new ClusterSumExternalField (new ClusterBonds[]{clusterbonds},new double[]{1},new MayerFunction[]{fRef});
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = new ClusterSum (new ClusterBonds[]{clusterbonds},new double[]{1},new MayerFunction[]{fRef});*/    
             
        /*int[][][] bondList1 =new int [][][] {{{0,1},{0,2}}};    	
        int[][][] bondList2 =new int [][][] {{{0,1},{1,2}}};
        int[][][] bondList3 =new int [][][] {{{0,1},{0,2},{1,2}}};    	  
        ClusterBonds cluster1 = new ClusterBonds(3, bondList1, false);
        ClusterBonds cluster2 = new ClusterBonds(3, bondList2, false);
        ClusterBonds cluster3 = new ClusterBonds(3, bondList3, false);    	
        double[] weights = {1.0/6.0, 1.0/3.0, 1.0/6.0};    	
        
        ClusterSumExternalField targetCluster = new ClusterSumExternalField (new ClusterBonds[] {cluster1, cluster2, cluster3},weights,new MayerFunction[]{fRef});
        ClusterAbstract refCluster = new ClusterSum (new ClusterBonds[]{cluster1, cluster2, cluster3},weights,new MayerFunction[]{fRef});*/
        
        /*int[][][] bondList1=new int [][][] {{{0,1},{1,2},{2,3}}};
        int[][][] bondList2=new int [][][] {{{0,1},{0,2},{2,3}}};
        int[][][] bondList3=new int [][][] {{{0,1},{1,2},{1,3}}};
        int[][][] bondList4=new int [][][] {{{0,1},{0,2},{0,3}}};
        int[][][] bondList5=new int [][][] {{{0,1},{0,2},{1,2},{1,3}}};
        int[][][] bondList6=new int [][][] {{{0,1},{0,2},{0,3},{1,2}}};
        int[][][] bondList7=new int [][][] {{{0,1},{1,2},{1,3},{2,3}}};
        int[][][] bondList8=new int [][][] {{{0,1},{1,2},{2,3},{0,3}}};
        int[][][] bondList9=new int [][][] {{{0,1},{0,2},{0,3},{1,2},{2,3}}};
        int[][][] bondList10=new int [][][] {{{0,1},{0,2},{1,2},{1,3},{2,3}}};
        int[][][] bondList11=new int [][][] {{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}};
        ClusterBonds cluster1 = new ClusterBonds(4, bondList1, false);
        ClusterBonds cluster2 = new ClusterBonds(4, bondList2, false);
        ClusterBonds cluster3 = new ClusterBonds(4, bondList3, false); 
        ClusterBonds cluster4 = new ClusterBonds(4, bondList4, false); 
        ClusterBonds cluster5 = new ClusterBonds(4, bondList5, false); 
        ClusterBonds cluster6 = new ClusterBonds(4, bondList6, false); 
        ClusterBonds cluster7 = new ClusterBonds(4, bondList7, false); 
        ClusterBonds cluster8 = new ClusterBonds(4, bondList8, false); 
        ClusterBonds cluster9 = new ClusterBonds(4, bondList9, false); 
        ClusterBonds cluster10 = new ClusterBonds(4, bondList10, false); 
        ClusterBonds cluster11 = new ClusterBonds(4, bondList11, false); 
        double[] weights = {1.0/4.0, 1.0/4.0, 1.0/8.0, 1.0/24.0, 1.0/4.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/24.0};    	
        ClusterSumExternalField targetCluster = new ClusterSumExternalField (new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8, cluster9, cluster10, cluster11},weights,new MayerFunction[]{fRef});
        ClusterAbstract refCluster = new ClusterSum (new ClusterBonds[]{cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8, cluster9, cluster10, cluster11},weights,new MayerFunction[]{fRef});*/
              
        ExternalVirialDiagrams refDiagrams = new ExternalVirialDiagrams(nPoints, false, false);
        ClusterSum refCluster = refDiagrams.makeRhoCluster(fRef, false);

        ExternalVirialDiagrams targetDiagrams = new ExternalVirialDiagrams(nPoints, false, false);
        ClusterSumExternalField targetCluster = (ClusterSumExternalField)targetDiagrams.makeRhoCluster(fRef, true);
        
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps");
        
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        final SimulationVirial sim = new SimulationVirial(space,new SpeciesFactorySpheres(), temperature,sampleCluster, refCluster,new ClusterAbstract[]{targetCluster});
        sim.setAccumulator(new AccumulatorRatioAverageCovariance());
        sim.box.getBoundary().setBoxSize(space.makeVector(new double []{10, 10, 10}));
        for (int i=0;i<nPoints; i++){
        	sim.box.getLeafList().getAtom(i).getPosition().setX(2, z0);
        }
       
        PotentialMaster potentialMaster=new PotentialMasterMonatomic (sim);
        P1HardBoundary p1wall=new P1HardBoundary (space);
        p1wall.setCollisionRadius(sigmaHSRef/2);
        p1wall.setActive(0, true, false);
        p1wall.setActive(0, false, false);
        p1wall.setActive(1, true, false);
        p1wall.setActive(1, false, false);
        potentialMaster.addPotential(p1wall, new IAtomType[]{sim.species.getAtomType(0)});
        targetCluster.setPotentialMaster(potentialMaster);
 
        sim.equilibrate(steps/40);
        
        System.out.println("equilibration finished");



        sim.ai.setMaxSteps(steps);
        
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
        
        sim.getController().actionPerformed();

        
       
      
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        System.out.println("ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("reference average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.ERROR.index)).getData()[0]);
        System.out.println("target average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverageCovariance.StatType.ERROR.index)).getData()[1]);

       
	}

    /**
     * Inner class for parameters
     */
    public static class VirialExternalFieldParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 1.5;
        public long numSteps = 100000;
        public double z0=-4.5;        
    }
}
