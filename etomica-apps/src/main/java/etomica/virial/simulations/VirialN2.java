/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.*;
import etomica.potential.P2SemiclassicalAtomic.AtomInfo;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

public class VirialN2 {
    public static void main(String[] args) {
        VirialN2Param params = new VirialN2Param();
        boolean isCommandLine = args.length > 0;
        if (isCommandLine) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
            // default options - choose these before committing to CVS
            //            params.potentialLevel = level.classical;
            //            params.temperatureK = 500;
            //            params.numSteps = (long)1E6;

            // runtime options - make changes in these and not the default options above
            params.nPoints = 3;
            params.potentialLevel = level.semiClassicalQFH;
            params.temperatureK = 500;
            params.numSteps = (long)1E6;
            params.nonAdditive = false;
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperatureK;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final boolean pairOnly = params.pairOnly;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        final boolean p2N2HellmannA = params.p2N2HellmannA;
        final boolean p3N2HellmannA = params.p3N2HellmannA;
        final boolean nonAdditive = params.nonAdditive;
        final level pLevel = params.potentialLevel;
        final int nT = params.nTimes;
        final boolean printInConvertedUnits = params.printInConvertedUnits;
        final double[] HSB = new double[8];

        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);

        System.out.println("Overlap sampling for N2 potential of Hellmann (2013) at " + temperatureK + " K");
        if (pLevel == level.semiClassicalQFH) System.out.println("Quadratic Feymann-Hibbs effective potential employed.");
        if (pLevel == level.semiClassicalTI) System.out.println("Takahashi-Imada effective potential employed.");
        System.out.println("Non - additive = "+nonAdditive);

        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");

        final Space space = Space3D.getInstance();
        // make ref and tar clusters
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = null;
        refCluster = new ClusterWheatleyHS(nPoints, fRef);

        ClusterAbstract tarCluster = null;

        final P2NitrogenHellmann p2N2 = new P2NitrogenHellmann(space);

        if (p2N2HellmannA) p2N2.parametersB = false;
        final P2SemiclassicalAtomic p2SCTI = new P2SemiclassicalAtomic(space, p2N2, temperature);
        final IPotentialAtomic p2 = (pLevel == level.classical ? p2N2 : (pLevel == level.semiClassicalQFH ? p2N2.makeQFH(temperature) : p2SCTI));

        PotentialMolecularMonatomic p2N2Molecular = new PotentialMolecularMonatomic(space, p2);
        MayerGeneral f2Tar = new MayerGeneral(p2N2Molecular);
        tarCluster = new ClusterWheatleySoft(nPoints, f2Tar, 1e-12);

        final double refIntegralF = HSB[nPoints];
        System.out.println("nPoints: "+nPoints);
        //        System.out.println("Ref Integral: HSB[2]");

        if (nPoints == 3 && nonAdditive) {
            final P3NitrogenHellmannNonAdditive p3N2NonAdditive = new P3NitrogenHellmannNonAdditive(space);
            if (p3N2HellmannA) p3N2NonAdditive.parametersB = false;
            final IPotentialAtomic p3 = p3N2NonAdditive;
            PotentialMolecularMonatomic p3N2Molecular = new PotentialMolecularMonatomic(space, p3);
            MayerFunctionMolecularThreeBody f3Tar = new MayerFunctionMolecularThreeBody(p3N2Molecular);
            MayerFunctionNonAdditive [] m1 = new MayerFunctionNonAdditive[]{null,null,null,f3Tar};
            tarCluster = new ClusterWheatleyMultibody(nPoints, f2Tar, m1);
        }
        refCluster.setTemperature(temperature);
        tarCluster.setTemperature(temperature);

        // make species
        ElementSimple n2 = new ElementSimple("N2", 2*Nitrogen.INSTANCE.getMass());
        final SpeciesSpheresRotating speciesN2 = new SpeciesSpheresRotating(space,n2);

        // make simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesN2, temperature, refCluster, tarCluster);
        // sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
        final Vector[] rv = space.makeVectorArray(4);
        rv[0].setX(0, massN*blN2*blN2*0.25);
        rv[0].setX(1, massN*blN2*blN2*0.25);
        p2SCTI.setAtomInfo(speciesN2.getAtomType(0), new AtomInfo() {
            @Override
            public Vector[] getMomentAndAxes(IAtomOriented molecule) {

                // rv[0,2] = 0
                // rv[3] is the orientation
                rv[3].E(molecule.getOrientation().getDirection());
                // rv[1] is an axis perpendicular to rv[3]
                rv[1].E(0);
                if (Math.abs(rv[3].getX(0)) < 0.5) {
                    rv[1].setX(0, 1);
                }
                else if (Math.abs(rv[3].getX(1)) < 0.5) {
                    rv[1].setX(1, 1);
                }
                else {
                    rv[1].setX(2, 1);
                }
                rv[2].Ea1Tv1(rv[1].dot(rv[3]), rv[3]);
                rv[1].ME(rv[2]);
                // rv[2] is an axis perpendicular to rv[3] and rv[1]
                rv[2].E(rv[1]);
                rv[2].XE(rv[3]);
                return rv;
            }
        });


        // add additional moves here, simulation already has translation and rotation moves

        System.out.println();
        String refFileName = params.refPrefFileName;
        //        if (isCommandLine) {
        //            String tempString = ""+temperatureK;
        //            if (temperatureK == (int)temperatureK) {
        //                // temperature is an integer, use "200" instead of "200.0"
        //                tempString = ""+(int)temperatureK;
        //            }
        //            refFileName = "refpref"+nPoints;
        //            refFileName += pairOnly ? "_2b" : "_3b";
        //            refFileName += "_"+tempString+"_";
        //            if (pLevel == level.semiClassical) refFileName += "SC";
        //            if (pLevel == level.classical) refFileName += "C";
        //        }

        long t1 = System.currentTimeMillis();
        // if using 3-body potential for B3, we must select initial configuration
        // that is not overlapping for any two molecules
        if (nPoints == 3 && nonAdditive) {
            IAtomList tarList = sim.box[1].getLeafList();
            for (int i = 0; i<tarList.size(); i++) {
                Vector p = tarList.get(i).getPosition();
                p.setX(i, 4.0);
            }
            sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
            sim.box[1].getSampleCluster().value(sim.box[1]);
        }

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
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        // this is where the simulation takes place

sim.getController().runActivityBlocking(new etomica.action.activity.ActivityIntegrate2(sim.integratorOS), 1000);
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
        double conv = 1.0;
        if (printInConvertedUnits) conv = Math.pow(1E-24/(Mole.UNIT.fromSim(1)),nPoints-1.0);

        double bn = ratio*refIntegralF*conv;
        double bnError = error*Math.abs(refIntegralF)*conv;

        System.out.println("ratio average: "+ratio+" error: "+error);
        System.out.println("abs average: "+bn+" error: "+bnError);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData();
        IData ratioData = allYourBase.getData(sim.accumulators[0].RATIO.index);
        IData ratioErrorData = allYourBase.getData(sim.accumulators[0].RATIO_ERROR.index);
        IData averageData = allYourBase.getData(sim.accumulators[0].AVERAGE.index);
        IData stdevData = allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(sim.accumulators[0].ERROR.index);
        IData correlationData = allYourBase.getData(sim.accumulators[0].BLOCK_CORRELATION.index);
        IData covarianceData = allYourBase.getData(sim.accumulators[0].BLOCK_COVARIANCE.index);
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
        ratioData = allYourBase.getData(sim.accumulators[1].RATIO.index);
        ratioErrorData = allYourBase.getData(sim.accumulators[1].RATIO_ERROR.index);
        averageData = allYourBase.getData(sim.accumulators[1].AVERAGE.index);
        stdevData = allYourBase.getData(sim.accumulators[1].STANDARD_DEVIATION.index);
        errorData = allYourBase.getData(sim.accumulators[1].ERROR.index);
        correlationData = allYourBase.getData(sim.accumulators[1].BLOCK_CORRELATION.index);
        covarianceData = allYourBase.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
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
    public static ClusterAbstract dummyTarCluster = null;
    public static ClusterAbstract dummyRefCluster = null;
    public static final double massN = Nitrogen.INSTANCE.getMass();
    public static final double blN2 = 1.1014;
    enum level {
        classical,semiClassicalQFH,semiClassicalTI;
    }
    /**
     * Inner class for parameters
     */
    public static class VirialN2Param extends ParameterBase {
        public int nPoints = 2;
        public double temperatureK = 500.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;
        public double sigmaHSRef = 4.50; // -1 means use equation for sigmaHSRef
        public level potentialLevel = level.classical;
        public boolean pairOnly = true;
        public boolean p2N2HellmannA = false;
        public boolean p3N2HellmannA = false;
        public boolean nonAdditive = false;
        public String jsonOutputFileName = "";
        public String refPrefFileName = null;
        public int nTimes = 1; // to run the simulation more than once
        public boolean printInConvertedUnits = false; // only prints bn and bn_error in converted units
    }
}
