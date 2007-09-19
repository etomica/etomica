package etomica.virial.simulations;



import etomica.action.Action;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.models.water.PNWaterGCPM;
import etomica.potential.Potential;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterCoupledFlipped;
import etomica.virial.ClusterSumPolarizable;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MCMoveClusterMoleculePushMulti;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryWaterGCPM;
import etomica.virial.cluster.Standard;


public class VirialWaterOverlapClusterCoupled extends Simulation {

    public static void main(String[] args) {

        int nPoints = 2;
        double temperature = Kelvin.UNIT.toSim(350);
        long steps = 10000l;
        int numSubSteps = 1000;
        double deltaDCut = 100;
        double pushR = 0;

        if (args.length > 0) nPoints = Integer.parseInt(args[0]);
        if (args.length > 1) temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[1]));
        if (args.length > 2) steps = Long.parseLong(args[2]);
        if (args.length > 3) deltaDCut = Double.parseDouble(args[3]);
        if (args.length > 4) pushR = Double.parseDouble(args[4]);

        double sigmaHSRef = 3.2;
        if (pushR > 2) {
            sigmaHSRef = pushR + 2;
        }
        final double[] HSB = new double[7];
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
//        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);

//        P2WaterSPCE pTarget = new P2WaterSPCE(space);
	//P2WaterTIP4P pTarget = new P2WaterTIP4P(space);
//        PotentialWaterPPC2 pTarget = new PotentialWaterPPC2(space);
//        PotentialWaterPPC9forB3 pTarget = new PotentialWaterPPC9forB3(space);
        final Potential pTarget = new PNWaterGCPM(space);
        
        // kmb added the code below; 9/23/05
        ClusterWeight sampleCluster1 = null;
     
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        ClusterAbstract targetCluster = Standard.virialClusterPolarizable(nPoints, fTarget, nPoints>3, false);
        ((ClusterSumPolarizable)targetCluster).setDeltaDCut(deltaDCut);
        targetCluster = new ClusterCoupledFlipped(targetCluster);
        if (targetCluster instanceof ClusterCoupledFlipped) {
            System.out.println("We're flipping");
        }

   	    sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster);

        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster);

       
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refSample.setTemperature(temperature);


        System.out.println(steps+" steps of size "+numSubSteps);

            final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactoryWaterGCPM(), temperature, new ClusterAbstract[]{refCluster,targetCluster},new ClusterWeight[]{refSample,sampleCluster1});


//            if (pushR > 0) {
                //((ClusterSumPolarizable)((ClusterCoupledFlipped)((ClusterWeightAbs)sim.meters[0].getClusters()[1]).getSubCluster()).getSubCluster()).pushR2 = pushR*pushR;
//                System.out.println("pushing to "+pushR);
//                sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveTranslate[1]);
//                MCMoveClusterMoleculePushMulti translateMove = new MCMoveClusterMoleculePushMulti(sim.integrators[1].getPotential(), sim.getRandom(), 1.0, nPoints-1);
//                translateMove.setMinRange(pushR);
//                sim.mcMoveTranslate[1] = translateMove;
//                sim.mcMoveTranslate[1] = new MCMoveClusterPullMulti(translateMove);
//                ((MCMoveClusterPullMulti)sim.mcMoveTranslate[1]).setMaxRange(200);
//                sim.integrators[1].getMoveManager().addMCMove(sim.mcMoveTranslate[1]);
                
//                sim.integratorOS.setAdjustStepFreq(false);
//                sim.integratorOS.setStepFreq0(0);
                
//                sim.initRefPref(null,steps/40);
//                pushR = 500;
//                System.out.println("pushing to "+pushR);
//                ((MCMoveClusterMoleculePushMulti)sim.mcMoveTranslate[1]).setMinRange(pushR);
//                
//                sim.initRefPref(null,steps/40);
//                pushR = 900;
//                System.out.println("pushing to "+pushR);
//                ((MCMoveClusterMoleculePushMulti)sim.mcMoveTranslate[1]).setMinRange(pushR);
//            }
            
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize((int)(steps/10));
            
//            final MeterDFVirial meterRDF = new MeterDFVirial(sim.getSpace()); //, sim.getSpeciesManager().getSpecies()[0]);
//            meterRDF.setBox(sim.box[1]);
//            meterRDF.getXDataSource().setXMax(1000);
//            meterRDF.getXDataSource().setNValues(500);
            
            
            for (int i=0; i<2; i++) {

                sim.integrators[i].getMoveManager().setEquilibrating(true);

                sim.mcMoveTranslate[i].setStepSize(10);

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
            String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+Kelvin.UNIT.fromSim(temperature) : null;
            sim.initRefPref(refFileName,steps/40);
            sim.equilibrate(refFileName,steps/20);
            if (pushR == 0 && (Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            //((ClusterSumPolarizable)((ClusterCoupledFlipped)sim.meters[1].getClusters()[0]).getSubCluster()).pushR2 = pushR*pushR;
//            sim.integrators[1].addIntervalAction(meterRDF);
//            sim.integrators[1].setActionInterval(meterRDF, 5);
            sim.setAccumulatorBlockSize((int)(steps*10));

            final int n = nPoints;
//            final ClusterAbstract tc = targetCluster;
//            Action rdfWriteAction = new Action() {
//                public void actionPerformed() {
//                    try {
//                        BoxCluster box = sim.box[1];
//                        CoordinatePairSet cPairs = box.getCPairSet();
//                        //System.out.println("writing xyz "+i);
//                        for (int j=0; j<n-1; j++) {
//                            for (int k=j+1; k<n; k++) {
//                                System.out.print(Math.sqrt(cPairs.getr2(j,k))+" ");
//                            }
//                        }
//                        System.out.println();
//                        DataGroup data = (DataGroup)meterRDF.getData();
//                        Data dataDF = data.getData(0);
//                        Data dataV = data.getData(1);
//                        Data dataV2 = data.getData(2);
//                        Data dataPi = data.getData(3);
//                        double deltaX = meterRDF.getXDataSource().getXMax() / meterRDF.getXDataSource().getNValues();
//                        int max = -1;
//                        for (int i=dataDF.getLength()-1; i>-1; i--) {
//                            if (dataDF.getValue(i) > 0) {
//                                max = i;
//                                break;
//                            }
//                        }
//                        if (max == -1) {
//                            System.out.println("oops nothing there!");
//                            return;
//                        }
//                        System.out.println("writing rdf");
//                        FileWriter fileWriter = new FileWriter("rdf"+n+(tc instanceof ClusterCoupledFlipped ? "f" : "")+(
//                                (pTarget instanceof PNWaterGCPM) ? "" : "_old")+".dat");
//                        for (int i=0; i<max+1; i++) {
//                            if (dataDF.getValue(i) > 0) {
//                                fileWriter.write((i+0.5)*deltaX+" "+dataDF.getValue(i)+" "+dataV.getValue(i)+" "+dataV2.getValue(i)+" "+dataPi.getValue(i)+"\n");
//                            }
//                        }
//                        fileWriter.close();
//                    }
//                    catch (IOException e) {throw new RuntimeException(e);}
//                }
//            };
            Action progressReport = new Action() {
                public void actionPerformed() {
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double ratio = sim.dsvo.getDataAsScalar();
                    double error = sim.dsvo.getError();
                    System.out.println("abs average: "+ratio*HSB[n]+", error: "+error*HSB[n]);
                }
            };
//            sim.integrators[1].addIntervalAction(rdfWriteAction);
//            sim.integrators[1].setActionInterval(rdfWriteAction, (int)(steps*10));
            sim.integratorOS.addIntervalAction(progressReport);
            sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

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

            //rdfWriteAction.actionPerformed();
            
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



