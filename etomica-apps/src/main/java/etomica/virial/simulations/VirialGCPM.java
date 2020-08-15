/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;


public class VirialGCPM {

    public static void main(String[] args) {
        VirialGCPMParam params = new VirialGCPMParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nPoints = 2;
            params.numSteps = 10000000;
            params.temperature = 300;
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        boolean additive = params.additive;
        int numSubSteps = 1000;
        double deltaCut = 100;

        double sigmaHSRef = 5;
        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Water overlap sampling B"+nPoints+" at T="+temperature);
        temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");

        Space space = Space3D.getInstance();

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        final IPotentialMolecular pTarget = new PNWaterGCPM(space);
        
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        ClusterAbstract targetCluster;
        if (additive) {
            // only additive (pairwise) contributions
            targetCluster = new ClusterWheatleySoft(nPoints, fTarget, 1e-12);
            if (nPoints == 2) {
                ((ClusterWheatleySoft)targetCluster).setDoCaching(false);
                targetCluster = new ClusterCoupledFlipped(targetCluster, space, 20);
            }
        }
        else {
            // full virial coefficient
            targetCluster = Standard.virialClusterPolarizable(nPoints, fTarget, nPoints>3, false);
            ((ClusterSumPolarizable)targetCluster).setDeltaCut(deltaCut);
            ((ClusterSumPolarizable)targetCluster).setCaching(false);
            targetCluster = new ClusterCoupledFlipped(targetCluster, space, 20);
        }

        ClusterWeight sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster);

        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster);

       
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refSample.setTemperature(temperature);

        SpeciesWater4P species = new SpeciesWater4P(space);
        species.setConformation(new ConformationWaterGCPM(space));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,
                temperature, new ClusterAbstract[]{refCluster,targetCluster},new ClusterWeight[]{refSample,sampleCluster1}, false);
        
        if  (targetCluster.value(sim.box[1]) == 0) {
            throw new RuntimeException("oops");
        }
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController2().addActivity(new ActivityIntegrate2(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }

        long t1 = System.currentTimeMillis();

        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
        sim.setAccumulatorBlockSize(steps);

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.integratorOS.setNumSubSteps((int)steps);
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize());

//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), 1000);
        long t2 = System.currentTimeMillis();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(HSB[nPoints]);
        System.out.println("time: "+(t2-t1)*0.001);
	}



    /**
     * Inner class for parameters
     */
    public static class VirialGCPMParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 350;   // Kelvin
        public long numSteps = 10000000;
        public boolean additive = true;
    }
}
