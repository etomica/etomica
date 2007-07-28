package etomica.virial.simulations;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.models.water.PotentialWaterGCPM3forB5;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterCoupled;
import etomica.virial.ClusterSumPolarizable;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryWater4P;
import etomica.virial.cluster.Standard;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;



/**

 * Generic simulation using Mayer sampling to evaluate cluster integrals

 */

public class VirialWaterOverlapClusterCoupled extends Simulation {





    public static void main(String[] args) {

        int nPoints = 2;

        double temperature = Kelvin.UNIT.toSim(350);

        long steps = 1000l;

        if (args.length > 0) nPoints = Integer.parseInt(args[0]);

        if (args.length > 1) temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[1]));

        if (args.length > 2) steps = Long.parseLong(args[2]);

        double sigmaHSRef = 3.2;

        double[] HSB = new double[7];

        HSB[2] = Standard.B2HS(sigmaHSRef);

        HSB[3] = Standard.B3HS(sigmaHSRef);

        HSB[4] = Standard.B4HS(sigmaHSRef);

        HSB[5] = Standard.B5HS(sigmaHSRef);

        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);

        System.out.println("B2HS: "+HSB[2]);

        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");

        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");

        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");

        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");

        System.out.println("Water Direct sampling B"+nPoints+" at T="+Kelvin.UNIT.fromSim(temperature));

		

        Space space = Space3D.getInstance();

        

        MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);

        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);

//        P2WaterSPCE pTarget = new P2WaterSPCE(space);
//	P2WaterTIP4P pTarget = new P2WaterTIP4P(space);
//        PotentialWaterPPC2 pTarget = new PotentialWaterPPC2(space);
//        PotentialWaterPPC9forB3 pTarget = new PotentialWaterPPC9forB3(space);
        //PotentialWaterGCPMforB3 pTarget = new PotentialWaterGCPMforB3(space);
        PotentialWaterGCPM3forB5 pTarget = new PotentialWaterGCPM3forB5(space);
  
        
        // kmb added the code below; 9/23/05
        ClusterWeight sampleCluster1 = null;
     
	    MayerGeneral fTarget = new MayerGeneral(pTarget);
	    MayerEGeneral eTarget = new MayerEGeneral(pTarget);
//	    ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, temperature);
	    ClusterSumPolarizable targetCluster = Standard.virialClusterPolarizable(nPoints, fTarget, nPoints>3, eTarget, false);
	    //	ClusterCoupled targetClusterCoupled = new ClusterCoupled(targetCluster);

// old "trunc" code before flipping molecules; KMB and AJS, 7/25/07
/*	    if (args.length > 3 && args[3].equals("trunc")) {
            Potential2Transformed pSample = new Potential2Transformed(space,pTarget);
            pSample.setAtomPositionDefinition(new AtomPositionFirstAtom());
            MayerGeneral fSample = new MayerGeneral(pSample);
            MayerEGeneral eSample = new MayerEGeneral(pSample);
            ClusterAbstract sampleCluster2 = Standard.virialCluster(nPoints, fSample, true, eSample, temperature); //, true);
            sampleCluster1 = ClusterWeightAbs.makeWeightCluster(sampleCluster2);
        }
        else {
*/        	    sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster.makeCopy());
  //      	    sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetClusterCoupled.makeCopy());
                             
//        }

        

        


        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, false, false);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster.makeCopy());

       
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refCluster.setTemperature(temperature);


//        System.out.println(steps+" steps of size "+defaults.blockSize);

		

//		while (true) {

            SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactoryWater4P(), temperature, new ClusterAbstract[]{refCluster,targetCluster},new ClusterWeight[]{refSample,sampleCluster1});
//            SimulationVirialOverlap sim = new SimulationVirialOverlap(space,defaults,new SpeciesFactoryWater(), temperature,refCluster,targetCluster);

/*            AtomTreeNodeWater secondWater = (AtomTreeNodeWaterGCPM)((Phase)sim.getPhaseList().get(1)).molecule(1).node;
            secondWater.H1.coord.position().PE(20);
            secondWater.H2.coord.position().PE(20);
            secondWater.M.coord.position().PE(20);
            secondWater.O.coord.position().PE(20);
            ((PhaseCluster)sim.getPhaseList().get(1)).trialNotify();
            ((PhaseCluster)sim.getPhaseList().get(1)).acceptNotify();
*/            
            
//            sim.integratorOS.setStepFreq0(0);
//            sim.integratorOS.setAdjustStepFreq(false);
            
            for (int i=0; i<2; i++) {

                sim.integrators[i].getMoveManager().setEquilibrating(true);

                sim.mcMoveTranslate[i].setStepSize(1.0);

//                sim.mcMoveTranslate[i].setAcceptanceTarget(0.20);

//                sim.mcMoveTranslate[i].setAdjustInterval(sim.mcMoveTranslate[i].getAdjustInterval()*100);

                sim.mcMoveRotate[i].setStepSize(1.5);

  //              sim.mcMoveRotate[i].setAdjustInterval(sim.mcMoveRotate[i].getAdjustInterval()*100);

//                sim.mcMoveRotate[i].setAcceptanceTarget(0.20);

               /* if (sim.mcMoveMulti != null) {

                    sim.mcMoveMulti[i].setStepSize(0.5);

                    sim.mcMoveMulti[i].setAdjustInterval(sim.mcMoveMulti[i].getAdjustInterval()*100);

                    sim.mcMoveMulti[i].setAcceptanceTarget(0.20);

                }*/

            }
            


//            sim.mcMoveAtom1.setNoisyAdjustment(true);

//            sim.mcMoveRotate.setNoisyAdjustment(true);

//            sim.mcMoveMulti.setNoisyAdjustment(true);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
            sim.initRefPref(refFileName,steps/100);
            sim.equilibrate(refFileName,steps/40);

/*            try { 

                fileReader = new FileReader("refpref");

                BufferedReader bufReader = new BufferedReader(fileReader);

                String refPrefString = bufReader.readLine();

                refPref = Double.parseDouble(refPrefString);

                bufReader.close();

                fileReader.close();

                System.out.println("setting ref pref to "+refPref);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,1),0);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,1),1);

                sim.setRefPref(refPref,1);

            }

            catch (IOException e) {

                // file not there, which is ok.

            }



            if (refPref == -1) {

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,101),0);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,101),1);

                sim.setRefPref(10000,15);

                sim.ai.setMaxSteps(steps/100);

                sim.getController().run();



                int newMinDiffLoc = sim.dsvo.minDiffLocation();

                int nBennetPoints = sim.accumulators[0].getNBennetPoints();

                refPref = sim.accumulators[0].getBennetBias(nBennetPoints-newMinDiffLoc-1);

                System.out.println("setting ref pref to "+refPref);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,11),0);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,11),1);

                sim.setRefPref(refPref,0.2);

                // reset to -1 so we know later on we're still searching for the right value

                refPref = -1;

                for (int i=0; i<2; i++) {

                    // call setPhase again so meter can re-synchronize with samplecluster

                    try {

                        sim.integrators[i].reset();

                    }

                    catch (ConfigurationOverlapException e) {}

                    sim.meters[i].setPhase(sim.phase[i]);

                }

                sim.getController().addAction(sim.ai);

            }



            // run a short simulation to get reasonable MC Move step sizes and

            // (if needed) narrow in on a reference preference

            sim.ai.setMaxSteps(steps/100);



            sim.getController().run();



            if (refPref == -1) {

                int newMinDiffLoc = sim.dsvo.minDiffLocation();

                int nBennetPoints = sim.accumulators[0].getNBennetPoints();

                refPref = sim.accumulators[0].getBennetBias(nBennetPoints-newMinDiffLoc-1);

                System.out.println("setting ref pref to "+refPref);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,1),0);

                sim.setAccumulator(new AccumulatorVirialOverlapSingleAverage(defaults.blockSize,1),1);

                sim.setRefPref(refPref,1);

                try {

                    FileWriter fileWriter = new FileWriter("refpref");

                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);

                    bufWriter.write(String.valueOf(refPref)+"\n");

                    bufWriter.close();

                    fileWriter.close();

                }

                catch (IOException e) {

                    throw new RuntimeException("couldn't write to refpref file");

                }

            }

            for (int i=0; i<2; i++) {

                sim.integrators[i].setEquilibrating(false);

                // call setPhase again so meter can re-synchronize with samplecluster

                try {

                    sim.integrators[i].reset();

                }

                catch (ConfigurationOverlapException e) {}

                sim.meters[i].setPhase(sim.phase[i]);

            }
*/
            sim.ai.setMaxSteps(steps);

            System.out.println("equilibration finished");

            for (int i=0; i<2; i++) {

                System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "

                                                        +sim.mcMoveRotate[i].getStepSize());

            }


            

            sim.getController().actionPerformed();



            System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());

            

            double ratio = sim.dsvo.getDataAsScalar();
            double error = sim.dsvo.getError();
            System.out.println("ratio average: "+ratio+", error: "+error);
            System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);

            DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
            System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
            System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
            System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

            allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
            System.out.println("water ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
            System.out.println("water average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
            System.out.println("water overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

            

/*		double[] valueArray = targetCluster.getMaxValueArray();
            
            	System.out.println("Here are the pi versus log O-O values");
            
            	for(int y=0; y<9; y++) {
            		System.out.println(valueArray[y]);
    			}
            
            	double[] binCountArray = targetCluster.getBinCountArray();
            
            	System.out.println("Here are the bin counts");
            
            	for(int y=0; y<10; y++) {
            		System.out.println(binCountArray[y]);
    			}
*/            

//       }

	}

}



