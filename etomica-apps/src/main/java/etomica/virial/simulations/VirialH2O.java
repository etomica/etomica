/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.data.IData;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.mcmove.MCMove;
import etomica.math.DoubleRange;
import etomica.models.water.P2WaterSzalewicz;
import etomica.models.water.P2WaterSzalewicz.Component;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialMolecular;
import etomica.potential.P2WaterPotentialsJankowski;
import etomica.potential.PotentialMolecularMonatomic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Kelvin;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import org.json.simple.JSONObject;

import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

//import etomica.virial.PotentialCommonAtomic;

public class VirialH2O {
    public static void main(String[] args) {
        VirialH2OParam params = new VirialH2OParam();
        boolean isCommandLine = args.length > 0;
        if (isCommandLine) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
            // default options - choose these before committing to CVS
            params.iSurf = 1;
            params.temperatureK = 500;
            params.numSteps = (long)1E6;

            // runtime options - make changes in these and not the default options above
//            params.nPoints = 3;
//            params.iSurf = 1;
//            params.temperatureK = 1000;
//            params.numSteps = (long)1E6;
//            params.nonAdditive = true;
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperatureK;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        final int iSurf = params.iSurf;        
        final int iMon = params.iMon;
        final boolean nonAdditive = params.nonAdditive;
        final double[] HSB = new double[8];

        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
//        double a = 2892.33657322022*1000*Calorie.UNIT.toSim(1)/Mole.UNIT.toSim(1);
//        System.out.println(a);
//        System.exit(1);
        System.out.println("iMon "+iMon+"; iSurf "+iSurf);        
        if (iSurf == 1) {
            System.out.println("Potential surface used : CCpol-8sf[2014] (Eckart embedding)");
        }
        else if (iSurf == 2) {
            System.out.println("Potential surface used : CCpol-8sfIR[2014] (Eckart embedding)");
        }
        else if (iSurf == 3) {
            System.out.println("Potential surface used : CCpol-8sfIR[2012] (Radau f=1 embedding)");
        }
        else if (iSurf == 4) {
            System.out.println("Potential surface used : CCpol-8sfIR[2012] (Eckart embedding)");
        }
        else if (iSurf == 5) {
            System.out.println("Potential surface used : SAPT-5s'f[2014]");
        }
        else if (iSurf == 6) {
            System.out.println("Potential surface used : SAPT-5s'fIR[2014]");
        }
        else if (iSurf == 7) {
            System.out.println("Potential surface used : SAPT-5s'fIR[2012]");
        }
        else if (iSurf == 8) {
            System.out.println("Potential surface used : SAPT-5s'f[2006]");
        }
        else if (iSurf == 9) {
            System.out.println("Potential surface used : SAPT-5s'fIR[2006]");
        }
        else {
            System.out.println("Potential surface used : CCpol-8sfIR[2014] (Radau f=1 embedding)");
        }
        System.out.println("Temperature : "+temperatureK+"K");
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
        
        final Space space = Space3D.getInstance();
        
        // make ref and tar clusters
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = new ClusterWheatleyHS(nPoints, fRef);

        final IPotentialAtomic p2H2ORigidSC = new P2WaterSzalewicz(space,2).makeSemiclassical(temperature);
//        final PotentialCommonAtomic p2H2ORSCCommon = new PotentialCommonAtomic(p2H2ORSC);
//        PotentialMolecularMonatomic p2H2ORSCMolecular = new PotentialMolecularMonatomic(space, p2H2ORSCCommon);
        PotentialMolecularMonatomic p2H2OMolecularRigidSC = new PotentialMolecularMonatomic(space, p2H2ORigidSC);
        
//        final IPotentialAtomic p2H2O = new P2WaterPotentialsJankowski(space,iSurf,iMon,temperature, p2H2ORSCCommon);
        final IPotentialAtomic p2H2OFlex = new P2WaterPotentialsJankowski(space,iSurf,iMon,temperature, p2H2ORigidSC);
        PotentialMolecularMonatomic p2H2OMolecularFlex = new PotentialMolecularMonatomic(space, p2H2OFlex);
        
        MayerGeneral f2TarFlex = new MayerGeneral(p2H2OMolecularFlex);        
        ClusterAbstract tarCluster = null;
        
        // difference cluster to compute difference between hybrid and rigid potentials
        MayerGeneral f2TarSubtract = new MayerGeneral(p2H2OMolecularRigidSC);
        ClusterAbstract tarSubtract = null;

        tarCluster = new ClusterWheatleySoft(nPoints, f2TarFlex, 1e-12);
        ((ClusterWheatleySoft)tarCluster).setDoCaching(false);

        tarSubtract = new ClusterWheatleySoft(nPoints, f2TarSubtract, 1e-12);
        ((ClusterWheatleySoft)tarSubtract).setDoCaching(false);
        
        if (nPoints == 3 && nonAdditive) {
            P2WaterSzalewicz p23cH2O = new P2WaterSzalewicz(space, 2);
            p23cH2O.setComponent(Component.NON_PAIR);            
            IPotentialAtomic p23cH2Oa = p23cH2O.makeSemiclassical(temperature);
//            PotentialCommonAtomic p23cH2OCommon = new PotentialCommonAtomic(p23cH2Oa);
            PotentialMolecularMonatomic p23H2O = new PotentialMolecularMonatomic(space, p23cH2Oa);
            P2WaterSzalewicz p3cH2O = new P2WaterSzalewicz(space, 3);
            p3cH2O.setComponent(Component.NON_PAIR);
            IPotentialAtomic p3aH2O = p3cH2O.makeSemiclassical(temperature);
//            PotentialCommonAtomic p3aH2OCommon = new PotentialCommonAtomic(p3aH2O);
            PotentialMolecularMonatomic p3H2O = new PotentialMolecularMonatomic(space, p3aH2O);
            MayerFunctionMolecularThreeBody f3H2O = new MayerFunctionMolecularThreeBody(new PotentialNonAdditive(new IPotentialMolecular[]{p23H2O,p3H2O}));
            MayerFunctionNonAdditive[] fNA = new MayerFunctionNonAdditive[] {null,null,null,f3H2O};
            
            tarCluster = new ClusterWheatleyMultibody(nPoints, f2TarFlex, fNA);
            ((ClusterWheatleyMultibody)tarCluster).setDoCaching(false);
            ((ClusterWheatleyMultibody)tarCluster).setRCut(100);
            
            tarSubtract = new ClusterWheatleyMultibody(nPoints, f2TarSubtract, fNA);
            ((ClusterWheatleyMultibody)tarSubtract).setDoCaching(false);
            ((ClusterWheatleyMultibody)tarSubtract).setRCut(100);
        }
        
        final ClusterDifference diffCluster = new ClusterDifference(tarCluster,new ClusterAbstract[] {tarSubtract});
        
        // pure B2 for water.  we need flipping.
        // additive B3 for water should be fine and biconnectivity will help with mixture coefficients.        
        double flipDist = nPoints==2? 20 : 10;
        final ClusterAbstract tarFlipped = new ClusterCoupledAtomFlipped(tarCluster,space,flipDist);
        final ClusterAbstract tarSubFlipped = new ClusterCoupledAtomFlipped(tarSubtract,space,flipDist);
        final ClusterAbstract diffClusterNew = new ClusterCoupledAtomFlipped(diffCluster, space, flipDist);
        System.out.println("Calulating difference from Cc-pol2-sc");
        System.out.println("Including flip moves starting at a distance of: "+flipDist);

        final double refIntegralF = HSB[nPoints];
        System.out.println("nPoints: "+nPoints+" Reference hard sphere diameter used: "+params.sigmaHSRef);
        if (nonAdditive) System.out.println("Non additive");
        else if (nPoints == 3) System.out.println("Additive");
        System.out.println("HSB["+nPoints+"]: "+refIntegralF);
        tarFlipped.setTemperature(temperature);
        tarSubFlipped.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        tarCluster.setTemperature(temperature);
        tarSubtract.setTemperature(temperature);
        diffClusterNew.setTemperature(temperature);
        // make species
        final SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O",Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        speciesH2O.setAxisSymmetric(false);
        
        // make simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesH2O, temperature, refCluster, diffClusterNew);
        int [] seeds = sim.getRandomSeeds();
        System.out.println("Random seeds: "+Arrays.toString(seeds));

        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
        

        // add additional moves here; simulation already has translation and rotation moves
        
        System.out.println();
        String refFileName = null;
        if (isCommandLine) {
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints;
            refFileName += (nPoints == 2) ? "_2b" : "_3b";
            refFileName += "_"+tempString+"_";
            refFileName += "C";            
        }       
        
        long t1 = System.currentTimeMillis(); // Start time for simulation
        
        // if using 3-body potential for B3, we must select initial configuration
        // that is not overlapping for any two molecules
        IAtomList tarList = sim.box[1].getLeafList();
        for (int i=0; i<tarList.getAtomCount(); i++) {
            Vector p = tarList.getAtom(i).getPosition();
            p.setX(i, 4.0);                
        }            
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
        sim.box[1].getSampleCluster().value(sim.box[1]);        
        sim.initRefPref(refFileName, steps/40);
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }
        
        sim.equilibrate(refFileName, steps/10);
        System.out.println("equilibration finished");        
        
        // Collecting histogram in reference system
        final HistogramNotSoSimple h1 = new HistogramNotSoSimple(new DoubleRange(0,100));        
        final HistogramNotSoSimple h2 = new HistogramNotSoSimple(new DoubleRange(0,100));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorInitialized(IntegratorEvent e) {}
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorStepFinished(IntegratorEvent e) {
                IAtomList atoms = sim.box[1].getLeafList();
                double x01 = Math.sqrt(atoms.getAtom(0).getPosition().Mv1Squared(atoms.getAtom(1).getPosition()));
                double x02 = Math.sqrt(atoms.getAtom(0).getPosition().Mv1Squared(atoms.getAtom(2).getPosition()));
                double x12 = Math.sqrt(atoms.getAtom(1).getPosition().Mv1Squared(atoms.getAtom(2).getPosition()));
                double xMax = Math.max(x01,x02);
                xMax = Math.max(xMax,x12);
                double y1 = tarFlipped.makeCopy().value(sim.box[1]);
                h1.addValue(xMax,Math.abs(y1));
                double y2 = tarSubFlipped.makeCopy().value(sim.box[1]);
                h2.addValue(xMax,Math.abs(y2));
            }
        };
        sim.integrators[1].getEventManager().addListener(histListenerTarget);
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        // this is where the simulation takes place
        sim.getController().actionPerformed();
        
        long t2 = System.currentTimeMillis(); // End time for simulation        
        
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
        
        double bn = ratio*refIntegralF;
        double bnError = error*Math.abs(refIntegralF);
        
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
        
        if (isCommandLine && params.jsonOutputFileName != null) {        
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
                jsonFile.write(JSONObject.toJSONString(resultsMap));
                jsonFile.write("\n");
                jsonFile.close();
            } catch (IOException e) {
                throw new RuntimeException(e.getMessage());
            }                
        }
        double [] x1 = h1.xValues();
        double [] hist1 = h1.getHistogram();
        for (int i=0; i<x1.length; i++) {
            if (!Double.isInfinite(hist1[i]) && !Double.isNaN(hist1[i])) {
                System.out.println(x1[i]+" "+hist1[i]);
            }
        }
        System.out.println("******");
        double [] x2 = h2.xValues();
        double [] hist2 = h2.getHistogram();
        for (int i=0; i<x2.length; i++) {
            if (!Double.isInfinite(hist2[i]) && !Double.isNaN(hist2[i])) {
                System.out.println(x2[i]+" "+hist2[i]);
            }
        }
    }
    /**
     * Inner class for parameters
     */
    public static class VirialH2OParam extends ParameterBase {
        public int nPoints = 2;        
        public double temperatureK = 500.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public double sigmaHSRef = 5; // -1 means use equation for sigmaHSRef        
        public int iSurf = 1, iMon = 0; // Parameters for the potential class
        public String jsonOutputFileName = null;
        public int nTimes = 1; // to run the simulation more than once
        public boolean nonAdditive = false;
    }
}
