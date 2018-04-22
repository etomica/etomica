/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.P2O2Bartolomei;
import etomica.potential.PotentialMolecularMonatomic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.ClusterWheatleySoft;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;

public class VirialO2 {
    public static void main(String[] args) {
        VirialO2Param params = new VirialO2Param();
        boolean isCommandLine = args.length > 0;
        if (isCommandLine) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
            // default options - choose these before committing to CVS
            params.nPoints = 2;
            params.potentialLevel = level.classical;
            params.s = 0;
            params.isPT2 = false;
            params.temperatureK = 500;
            params.numSteps = (long)1E6;            

            // runtime options - make changes in these and not the default options above
//            params.nPoints = 2;
//            params.potentialLevel = level.classical;
//            params.temperatureK = 600;
//            params.numSteps = (long)1E6;
//            params.s = 0;
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperatureK;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final boolean pairOnly = params.pairOnly;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;        
        final level pLevel = params.potentialLevel;
        final int nT = params.nTimes;
        final int s = params.s;
        final boolean isPT2 = params.isPT2;
        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        
        System.out.println("Overlap sampling for O2 pair potential of Bartolomei et al. (2010) at " + temperatureK + " K");        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
        System.out.println("multiplicity: "+s+" PT2: "+isPT2);
        
        Space space = Space3D.getInstance();
        // make ref and tar clusters
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);
        ClusterAbstract tarCluster = null;        
                
        P2O2Bartolomei p2 = new P2O2Bartolomei(space);
        P2O2Bartolomei.setS(s);
        P2O2Bartolomei.setPT2(isPT2);
        PotentialMolecularMonatomic p2O2Molecular = new PotentialMolecularMonatomic(space, p2);        
        MayerGeneral f2Tar = new MayerGeneral(p2O2Molecular);
        tarCluster = new ClusterWheatleySoft(nPoints, f2Tar, 1e-12);
        tarCluster.setTemperature(temperature);
        
        // make species
        SpeciesSpheresRotating speciesUranium = new SpeciesSpheresRotating(space,new ElementSimple("U",238.02891));
        
        // make simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesUranium, temperature, refCluster, tarCluster);
//        sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
        
        // add additional moves here, simulation already has translation and rotation moves
        
        System.out.println();
        String refFileName = null;
        if (isCommandLine) {
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints;
            refFileName += pairOnly ? "_2b" : "_3b";
            refFileName += "_"+tempString+"_";
            if (pLevel == level.classical) refFileName += "C";
        }
        long t1 = System.currentTimeMillis();
        // if using 3-body potential for B3, we must select initial configuration
        // that is not overlapping for any two molecules        
        
        sim.initRefPref(refFileName, steps/40);
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }
        
        sim.equilibrate(refFileName, steps/10);
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }
        final double refIntegralF = HSB[nPoints];
        if (!isCommandLine) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    //                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }
        // this is where the simulation takes place
        sim.getController().actionPerformed();
        //end of simulation
        long t2 = System.currentTimeMillis();
        
        
        System.out.println(" ");
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        System.out.println(" ");
        
        System.out.println("Reference system: ");
        List<MCMove> refMoves = sim.integrators[0].getMoveManager().getMCMoves();
        for (MCMove m : refMoves) {
            double acc = m.getTracker().acceptanceRatio();
            System.out.println(m.toString()+" acceptance ratio: "+acc);            
        }
        System.out.println("Target system: ");
        List<MCMove> tarMoves = sim.integrators[1].getMoveManager().getMCMoves();
        for (MCMove m : tarMoves) {
            double acc = m.getTracker().acceptanceRatio();

//            if (acc == 1) {
//                throw new RuntimeException("something seems fishy");
//            }
            System.out.println(m.toString()+" acceptance ratio: "+acc);
        }

        // Printing results here
        double[] ratioAndError = sim.dvo.getAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        double bn = ratio*HSB[nPoints];
        double bnError = error*Math.abs(HSB[nPoints]);
        
        System.out.println("ratio average: "+ratio+" error: "+error);
        System.out.println("abs average: "+bn+" error: "+bnError);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData();
        IData ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
        IData ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
        IData averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        IData stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        IData correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        IData covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(3));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        double refAvg = averageData.getValue(0);
        double refOvAvg = averageData.getValue(1);
        System.out.print(String.format("reference ratio average: %20.15e error:  %10.5e  cor: %6.4f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
        System.out.print(String.format("reference average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("reference overlap average: %20.15e stdev: %9.4e error: %9.3e cor: %6.4f\n",
                              averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));
        
        allYourBase = (DataGroup)sim.accumulators[1].getData();
        ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
        ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
        averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        int n = sim.numExtraTargetClusters;
        correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        double tarAvg = averageData.getValue(0);
        double tarOvAvg = averageData.getValue(1);
        double tarCorr = correlationData.getValue(0);
        System.out.print(String.format("target ratio average: %20.15e  error: %10.5e  cor: %6.4f\n", ratioData.getValue(n+1), ratioErrorData.getValue(n+1), correlationCoef));
        System.out.print(String.format("target average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("target overlap average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
                              averageData.getValue(n+1), stdevData.getValue(n+1), errorData.getValue(n+1), correlationData.getValue(n+1)));
        if (!isCommandLine) {
            if ((t2-t1)/1000.0 > 24*3600) {            
                System.out.println("time: "+(t2-t1)/(24*3600*1000.0)+" days");
            }
            else if ((t2-t1)/1000.0 > 3600) {            
                System.out.println("time: "+(t2-t1)/(3600*1000.0)+" hrs");
            }
            else if ((t2-t1)/1000.0 > 60) {           
                System.out.println("time: "+(t2-t1)/(60*1000.0)+" mins");
            }
            else {           
                System.out.println("time: "+(t2-t1)/1000.0+" secs");
            }
        }
        if (isCommandLine) {
            LinkedHashMap resultsMap = new LinkedHashMap();
            resultsMap.put("temperature", temperatureK);
            resultsMap.put("bn", bn);
            resultsMap.put("bnError", bnError);
            resultsMap.put("refAvg", refAvg);
            resultsMap.put("refOvAvg", refOvAvg);
            resultsMap.put("tarAvg", tarAvg);
            resultsMap.put("tarOvAvg", tarOvAvg);
            resultsMap.put("tarCorr", tarCorr);

            if ((t2-t1)/1000.0 > 24*3600) {
                resultsMap.put("time",(t2-t1)/(24*3600*1000.0));
                resultsMap.put("unit","days");
                System.out.println("time: "+(t2-t1)/(24*3600*1000.0)+" days");
            }
            else if ((t2-t1)/1000.0 > 3600) {
                resultsMap.put("time",(t2-t1)/(3600*1000.0));
                resultsMap.put("unit","hrs");
                System.out.println("time: "+(t2-t1)/(3600*1000.0)+" hrs");
            }
            else if ((t2-t1)/1000.0 > 60) {
                resultsMap.put("time",(t2-t1)/(60*1000.0));
                resultsMap.put("unit","mins");
                System.out.println("time: "+(t2-t1)/(60*1000.0)+" mins");
            }
            else {
                resultsMap.put("time",(t2-t1)/(1000.0));
                resultsMap.put("unit","secs");
                System.out.println("time: "+(t2-t1)/1000.0+" secs");
            }                   

            try {
                FileWriter jsonFile = new FileWriter(params.jsonOutputFileName);
                ObjectMapper om = new ObjectMapper();
                jsonFile.write(om.writeValueAsString(resultsMap));
                jsonFile.write("\n");
                jsonFile.close();
            } catch (IOException e) {
                throw new RuntimeException(e.getMessage());
            }                
        }

    }
    enum level {
        classical
    }
    /**
     * Inner class for parameters
     */
    public static class VirialO2Param extends ParameterBase {
        public int nPoints = 2;        
        public double temperatureK = 500.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public double sigmaHSRef = 4.50; // -1 means use equation for sigmaHSRef        
        public level potentialLevel = level.classical;
        public boolean pairOnly = true;        
        public String jsonOutputFileName = "";
        public int nTimes = 1; // to run the simulation more than once
        public int s = 0; // multiplicity
        public boolean isPT2 = false;
    }

}
