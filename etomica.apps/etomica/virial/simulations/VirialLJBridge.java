package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSumEF;
import etomica.virial.ConfigurationClusterMove;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * MSMC computations (via overlap sampling) of the bridge diagram(s) within the f-bond-only formulation of B4 and B5
 * 
 * Kate Shaul
 */
public class VirialLJBridge {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();

        double temperature; final int nPoints; double sigmaHSRef;
        long steps;
        if (args.length == 0) {
        	
        	nPoints = params.numMolecules;
            temperature = params.temperature;
            steps = params.numSteps;
            sigmaHSRef = params.sigmaHSRef;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 4) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperature = Double.parseDouble(args[1]);
            steps = Integer.parseInt(args[2]);
            sigmaHSRef = Double.parseDouble(args[3]);
            params.writeRefPref = true;
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed to VirialRowleyAlcohol.");
        }

        

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
        System.out.println("Lennard Jones overlap sampling f-bond decomposition of B"+nPoints+" at T= "+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        MayerESpherical eTarget = new MayerESpherical(pTarget);
        
        ClusterAbstract targetCluster;
        
        if (nPoints == 4) {
       
       
        
	        int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0},{0,2},{1,3}}, {}};
	        
        	
        	//ring
	        //int[][][] bondList = {{{0,1},{1,2},{2,3},{3,0}}, {}};
	        
	        ClusterBonds cluster = new ClusterBonds(4, bondList, false);
	        
	        
	
	        double[] weights = {-1.0/8.0};
	        
	        //double[] weights = {-3.0/8.0};
	  
		    
		    targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster},weights,new MayerFunction[]{eTarget});
        } else if (nPoints == 5) {
        	
        	int[][][] bondList1 = {{{0,1},{0,2},{0,4},{1,2},{1,4},{2,3},{3,4}}, {}};
        	
        	int[][][] bondList2 = {{{0,1},{0,2},{0,4},{1,2},{1,4},{2,3},{3,4},{2,4}}, {}};
        	
        	int[][][] bondList3 = {{{0,1},{0,2},{0,4},{1,2},{1,4},{2,3},{3,4},{1,3}}, {}};
        	
        	int[][][] bondList4 = {{{0,1},{0,2},{0,4},{1,2},{1,4},{2,3},{3,4},{1,3},{0,3}}, {}};
        	
        	int[][][] bondList5 = {{{0,1},{0,2},{0,4},{1,2},{1,4},{2,3},{3,4},{1,3},{0,3},{2,4}}, {}};
        	
        	ClusterBonds cluster1 = new ClusterBonds(5, bondList1, false);
        	ClusterBonds cluster2 = new ClusterBonds(5, bondList2, false);
        	ClusterBonds cluster3 = new ClusterBonds(5, bondList3, false);
        	ClusterBonds cluster4 = new ClusterBonds(5, bondList4, false);
        	ClusterBonds cluster5 = new ClusterBonds(5, bondList5, false);
        	
        	double[] weights = {-1.0, -1.0, -0.5, -1.0/3.0, -1.0/30.0};
        	
        	targetCluster =  new ClusterSumEF(new ClusterBonds[] {cluster1, cluster2, cluster3, cluster4, cluster5},weights,new MayerFunction[]{eTarget});
        
        
        } else {
        	throw new IllegalArgumentException("Cannot yet compute correction to f-bond decomposition approximation for that order of virial coefficient");
        }

        
        
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactorySpheres(), temperature,refCluster,targetCluster);
        
        
        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box[1]);
        
        sim.integratorOS.setNumSubSteps(1000);
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
                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
            }
        };
       // sim.integratorOS.addIntervalAction(progressReport);
       // sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("lennard jones ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("lennard jones average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("lennard jones overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
    	
    	// number of molecules in simulation (e.g., 2 for B2 calculation)

 
        public int numMolecules = 5;
        public double temperature = 1.0;
        public long numSteps = 10000;
        public double sigmaHSRef = 1.5;
        public boolean writeRefPref = false;
    }
    
}
