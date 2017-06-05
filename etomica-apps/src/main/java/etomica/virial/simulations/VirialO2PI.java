/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import etomica.integrator.IntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IMoleculeList;
import etomica.atom.AtomHydrogen;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Oxygen;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.P2O2Bartolomei;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.BohrRadius;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import org.json.simple.JSONObject;

import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

public class VirialO2PI {
    public static final double blO2 = BohrRadius.UNIT.toSim(2.28);
    public static final double massO = Oxygen.INSTANCE.getMass();

    public static void main(String[] args) {
        VirialO2Param params = new VirialO2Param();
        boolean isCommandLine = args.length > 0;
        if (isCommandLine) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
            // default options - choose these before committing to CVS
            params.nBeads = 8;
            params.temperatureK = 500;
            params.numSteps = (long)1E6;
            params.s = 0;
            params.isPT2 = false;

            // runtime options - make changes in these and not the default options above
//            params.nBeads = 8;
//            params.temperatureK = 1000;
//            params.s = 0;
//            params.isPT2 = true;
//            params.numSteps = (long)1E6;
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperatureK;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final boolean pairOnly = params.pairOnly;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        final int nBeads = params.nBeads;
        final int beadFac = params.beadFac;
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
        System.out.println("Path Integral Monte Carlo (PIMC) calculation with P = "+nBeads);
        System.out.println("multiplicity: "+s+" PT2: "+isPT2);


        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");

        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");

        Space space = Space3D.getInstance();
        // make ref and tar clusters
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        final P2O2Bartolomei p2 = new P2O2Bartolomei(space);
        P2O2Bartolomei.setS(s);
        P2O2Bartolomei.setPT2(isPT2);
        PotentialGroupPI pTarGroup = new PotentialGroupPI(beadFac);
        pTarGroup.addPotential(p2, new ApiIntergroupCoupled());
        MayerGeneral fTar = new MayerGeneral(pTarGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        ClusterWheatleySoft tarCluster = new ClusterWheatleySoft(nPoints, fTar, 1e-12);
        tarCluster.setTemperature(temperature);

        // make species
        AtomTypeOriented atype = new AtomTypeOriented(Oxygen.INSTANCE, space);
        SpeciesSpheresHetero speciesO2 = null;
        speciesO2 = new SpeciesSpheresHetero(space, new AtomTypeOriented[]{atype}) {
            protected IAtom makeLeafAtom(AtomType leafType) {
                return new AtomHydrogen(space, (AtomTypeOriented) leafType, blO2);
            }
        };
        speciesO2.setChildCount(new int [] {nBeads});
        speciesO2.setConformation(new ConformationLinear(space, 0));
        // make simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesO2, temperature, refCluster, tarCluster);
//        sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;

//        add additional moves here, simulation already has translation and rotation moves
//        rotation is a bit pointless when we can regrow the ring completely

        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
            sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        }

//        ring regrow translation
        MCMoveClusterRingRegrow refTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*massO*temperature);
        refTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

        MCMoveClusterRingRegrow tarTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        tarTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().addMCMove(refTr);
            sim.integrators[1].getMoveManager().addMCMove(tarTr);
        }
//        ring regrow orientation
        MCMoveClusterRingRegrowOrientation refOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
        MCMoveClusterRingRegrowOrientation tarOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
        refOr.setStiffness(temperature, massO);
        tarOr.setStiffness(temperature, massO);

        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().addMCMove(refOr);
            sim.integrators[1].getMoveManager().addMCMove(tarOr);
        }

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
            refFileName += "_"+tempString+"_PI";
        }
        long t1 = System.currentTimeMillis();

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

        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        final double refIntegralF = HSB[nPoints];
        if (! isCommandLine) {
            IIntegratorListener progressReport = new IIntegratorListener() {
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
            System.out.println(m.toString() + " acceptance ratio: " + acc);
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
        if (isCommandLine) {
            LinkedHashMap resultsMap = new LinkedHashMap();
            resultsMap.put("temperature", temperatureK);
            resultsMap.put("P", nBeads);
            resultsMap.put("bn", bn);
            resultsMap.put("bnError", bnError);
            resultsMap.put("refAvg", refAvg);
            resultsMap.put("refOvAvg", refOvAvg);
            resultsMap.put("tarAvg", tarAvg);
            resultsMap.put("tarOvAvg", tarOvAvg);
            resultsMap.put("tarCorr", tarCorr);
            for (MCMove m : tarMoves) {
                double acc = m.getTracker().acceptanceRatio();
                resultsMap.put(m.toString(), acc);
            }

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
//        sim.printResults(HSB[nPoints]);

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
    }
    
    /**
     * Inner class for parameters
     */
    public static class VirialO2Param extends ParameterBase {
        public int nPoints = 2;
        public int nBeads = 8;
        public double temperatureK = 500.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;        
        public double sigmaHSRef = 4.50; // -1 means use equation for sigmaHSRef
        public boolean pairOnly = true;        
        public int beadFac = 2;
        public int s = 0; // multiplicity
        public int nTimes = 1;
        public String jsonOutputFileName = "";
        public boolean isPT2 = false;
    }

}
