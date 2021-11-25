/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.Standard;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJMultiOverlap {


    public static void main(String[] args) {

        VirialMixParam params = new VirialMixParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // play with params here
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long numSteps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int mixID = params.mixID;
        int[] nTypes = params.nTypes;

        // sanity-check nPoints
        int sum = 0;
        for (int i=0; i<nTypes.length; i++) {
            sum += nTypes[i];
        }
        if (sum != nPoints) {
            throw new RuntimeException("Number of each type needs to add up to nPoints");
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
        System.out.println("Lennard Jones overlap sampling B"+nPoints+" at T="+temperature);
        System.out.println("Mixture type "+mixID);
        System.out.print("Points of each species:");
        for (int i=0; i<nTypes.length; i++) {
            System.out.print(" "+nTypes[i]);
        }
        System.out.println();
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        double sigma11 = 1.0;
        double sigma22 = 1.0;
        double sigma12 = 1.0;
        double epsilon11 = 1.0;
        double epsilon22 = 1.0;
        double epsilon12 = 1.0;
        if (mixID == 1) {
            epsilon12 = 0.75;
        }
        else if (mixID == 2) {
            sigma12 = 0.885;
            sigma22 = 0.769;
        }
        else if (mixID == 3) {
            sigma12 = 0.884; // who smoked crack?
            sigma22 = 0.768; // not me!
            epsilon12 = 0.773;
            epsilon22 = 0.597;
        }
        else if (mixID == 4) {
            epsilon12 = Math.sqrt(0.5);
            epsilon22 = 0.5;
        }
        else if (mixID != 0) {
            throw new RuntimeException("Don't know how to do mix "+mixID);
        }
        Potential2Soft p11Target = new P2LennardJones(sigma11, epsilon11);
        MayerGeneralSpherical f11Target = new MayerGeneralSpherical(p11Target);
        Potential2Soft p12Target = new P2LennardJones(sigma12, epsilon12);
        MayerGeneralSpherical f12Target = new MayerGeneralSpherical(p12Target);
        Potential2Soft p22Target = new P2LennardJones(sigma22, epsilon22);
        MayerGeneralSpherical f22Target = new MayerGeneralSpherical(p22Target);
        MayerESpherical e11Target = new MayerESpherical(p11Target);
        MayerESpherical e12Target = new MayerESpherical(p12Target);
        MayerESpherical e22Target = new MayerESpherical(p22Target);
        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{f11Target,f12Target},{f12Target,f22Target}},
                                                                               new MayerFunction[][]{{e11Target,e12Target},{e12Target,e22Target}}, nTypes);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((numSteps*1000)+" steps ("+numSteps+" blocks of 1000)");

        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")));
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints}, temperature, refCluster, targetCluster);
        sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(10*numSteps);
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, numSteps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, numSteps/40);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, numSteps);
        System.out.println("equilibration finished");

        if (false) {
            IAction progressReport = new IAction() {
                public void actionPerformed() {
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
            progressReportListener.setInterval((int)(numSteps/10));
            sim.integratorOS.getEventManager().addListener(progressReportListener);
        }

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        sim.getController().runActivityBlocking(ai);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
    }
    
    /**
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 5;
        public double temperature = 1.0;
        public long numSteps = 1000;
        public double sigmaHSRef = 1.5;
        public int mixID = 0;
        public int[] nTypes = new int[]{nPoints-1,1};
    }
}

