/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.models.water.*;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.Species;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialWaterOverlap {


    public static void main(String[] args) {
        int nPoints = 3;
        double temperature = Kelvin.UNIT.toSim(750);
        long steps = 50000l;
        int model = 0; //SPCE
        String[] modelStrings = new String[]{"SPCE","SPC","TIP4P"};
        String modelString = modelStrings[model];
        if (args.length > 0) nPoints = Integer.parseInt(args[0]);
        if (args.length > 1) temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[1]));
        if (args.length > 2) steps = Long.parseLong(args[2]);
        if (args.length > 3) {
            modelString = args[3];
            model = -1;
            for (int i=0; i<modelStrings.length; i++) {
                if (modelString.equals(modelStrings[i])) {
                    model = i;
                    break;
                }
            }
            if (model == -1) {
                throw new IllegalArgumentException("unknown model "+modelString);
            }
        }
        
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
        System.out.println("Water ("+modelString+") Overlap sampling B"+nPoints+" at T="+Kelvin.UNIT.fromSim(temperature));
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);

        IPotentialMolecular pTarget = null;
        ISpecies species = null;
        switch (model) {
            case 0: // SPCE
                pTarget = new P2WaterSPCE(space);
                species = SpeciesWater3P.create();
                break;
            case 1: // SPC
                pTarget = new P2WaterSPC(space);
                species = SpeciesWater3P.create();
                break;
            case 2: // TIP4P
                pTarget = new P2WaterTIP4P(space);
                species = SpeciesWater4P.create();
                break;
            default:
                throw new RuntimeException("unknown model "+model);
        }
        
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerEGeneral eTarget = new MayerEGeneral(pTarget);
        

        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, nPoints > 5);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, nPoints > 5);
        refCluster.setTemperature(temperature);
        
        ClusterWeight targetSampleCluster = ClusterWeightAbs.makeWeightCluster(targetCluster);
        ClusterWeight refSampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);

        if (nPoints == 2) {
            ((ClusterSum)targetCluster).setCaching(false);
            targetCluster = new ClusterCoupledFlipped(targetCluster, space);
            targetCluster.setTemperature(temperature);
            targetSampleCluster = new ClusterWeightAbs(targetCluster);
        }

        int blockSize = (int)steps/10;
        int numSubSteps = 1000;
        System.out.println(steps+" steps of size "+numSubSteps);
		
        SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{refSampleCluster, targetSampleCluster}, false);
        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(blockSize);
        
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, steps);
System.out.println("equilibration finished");

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+//sim.mcMoveAtom1[i].getStepSize()+" "
                                                    +sim.mcMoveRotate[i].getStepSize()
                                                    +((sim.mcMoveTranslate==null) ? "" : (" "+sim.mcMoveTranslate[i].getStepSize())));
        }
sim.getController().runActivityBlocking(ai);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        
        sim.printResults(HSB[nPoints]);
	}
}

