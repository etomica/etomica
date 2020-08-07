/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.models.OPLS.AceticAcidModPotentialHelper;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Mayer-sampling MC simulation for acetic acid using IMPROVED OPLS united-atom model
 *   
 *   
 *@author Hye Min Kim
 * Jan, 2011
 *  
 */

public class VirialAceticAcid {	
    public static void main(String[] args) {
    	
        VirialParam params = new VirialParam();
        Space space = Space3D.getInstance();
        ParseArgs.doParseArgs(params, args);
        
        final int nBody = params.nBody;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef=params.sigmaHSRef;
        
        final double[] HSB = new double[8];
        
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B5HS(sigmaHSRef);
        
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nBody+"HS: "+HSB[nBody]);
        System.out.println("IMPROVED OPLS Acetic acid B"+nBody+" at "+temperature+"K");     
        temperature = Kelvin.UNIT.toSim(temperature);//T * kB
		
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        
        PotentialGroup p = new PotentialGroup(2);

        MayerGeneral fTarget = new MayerGeneral(p);
        MayerEGeneral eTarget = new MayerEGeneral(p);
         
        ClusterAbstract targetCluster = Standard.virialCluster(nBody, fTarget, nBody>3, eTarget, true);
		if (nBody == 2){
			System.out.println("Flipping is applied");
			((ClusterSum)targetCluster).setCaching(false);
	        targetCluster = new ClusterCoupledFlipped(targetCluster, space);
		}
        targetCluster.setTemperature(temperature);
         
        ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody>3, eRef, true);
        refCluster.setTemperature(temperature);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;
        
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2 (space,new SpeciesAceticAcid(space),temperature,refCluster,targetCluster, false);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
    	SpeciesAceticAcid species = (SpeciesAceticAcid)sim.getSpecies(0);
    	AceticAcidModPotentialHelper.initPotential(space, species, p);

    	// bond angle bending is governed by harmonic potential
        double thetaEqCCDBO = 126*Math.PI/180;//equilibrium bond angle [=] degree, C-C=O
        double thetaEqDBOCSBO = 123*Math.PI/180;//O=C-O
        double thetaEqCSBOH = 107*Math.PI/180;//C-O-H
        double thetaEqCCSBO = 111*Math.PI/180;//C-C-O
        double kThetaCCDBO = Kelvin.UNIT.toSim(40300); //force constant [=] K, C-C=O
        double kThetaDBOCSBO = Kelvin.UNIT.toSim(40300);//O=C-O
        double kThetaCSBOH = Kelvin.UNIT.toSim(17600);//C-O-H
        double kThetaCCSBO = Kelvin.UNIT.toSim(35300);//C-C-O
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        
        P3BondAngle uBendingCCDBO = new P3BondAngle(space);//C-C=O
        P3BondAngle uBendingDBOCSBO = new P3BondAngle(space);//O=C-O
        P3BondAngle uBendingCSBOH = new P3BondAngle(space);//C-O-H
        P3BondAngle uBendingCCSBO = new P3BondAngle(space);//C-C-O
        
        uBendingCCDBO.setAngle(thetaEqCCDBO);
        uBendingCCDBO.setEpsilon(kThetaCCDBO);
        uBendingDBOCSBO.setAngle(thetaEqDBOCSBO);
        uBendingDBOCSBO.setEpsilon(kThetaDBOCSBO);
        uBendingCSBOH.setAngle(thetaEqCSBOH);
        uBendingCSBOH.setEpsilon(kThetaCSBOH);
        uBendingCCSBO.setAngle(thetaEqCCSBO);
        uBendingCCSBO.setEpsilon(kThetaCCSBO);
        
        pIntra.addPotential(uBendingCCDBO, new Atomset3IteratorIndexList(new int[][] {{0,1,2}}));//CH3:0, C:1, dBO:2, sBO:3, H:4
        pIntra.addPotential(uBendingDBOCSBO, new Atomset3IteratorIndexList(new int[][] {{2,1,3}}));
        pIntra.addPotential(uBendingCSBOH, new Atomset3IteratorIndexList(new int[][] {{1,3,4}}));
        pIntra.addPotential(uBendingCCSBO, new Atomset3IteratorIndexList(new int[][] {{0,1,3}}));
        sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{species});

        //utorsional = c1*[1+cos(phi+f1)]+c2*[1-cos(2*phi)], c1/kB = 630K, c2/kB = 1562.4K, f1=180°(for OCOH), 0°(for CCOH)

        MCMoveClusterTorsionAceticAcid[] torsionMoves1 = null;
        
        P4BondTorsion p4OCOH = new P4BondTorsion(space, 2*Kelvin.UNIT.toSim(630.0), Kelvin.UNIT.toSim(-630.0), Kelvin.UNIT.toSim(781.2), Kelvin.UNIT.toSim(0.0));
        P4BondTorsion p4CCOH = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(630.0), Kelvin.UNIT.toSim(781.2), Kelvin.UNIT.toSim(0.0));
        pIntra.addPotential(p4OCOH, new Atomset4IteratorIndexList(new int[][] {{2,1,3,4}}));
        pIntra.addPotential(p4CCOH, new Atomset4IteratorIndexList(new int[][] {{0,1,3,4}}));
        torsionMoves1 = new MCMoveClusterTorsionAceticAcid[2];
        torsionMoves1[0] = new MCMoveClusterTorsionAceticAcid(sim.integrators[1].getPotentialMaster(), space, sim.getRandom());
        sim.integrators[0].getMoveManager().addMCMove(torsionMoves1[0]);//reference system
        torsionMoves1[1] = new MCMoveClusterTorsionAceticAcid(sim.integrators[1].getPotentialMaster(), space, sim.getRandom());
        sim.integrators[1].getMoveManager().addMCMove(torsionMoves1[1]);//target system
        
        MCMoveClusterWiggleAceticAcid[] wiggle = new MCMoveClusterWiggleAceticAcid[2];
        wiggle[0] = new MCMoveClusterWiggleAceticAcid(sim,sim.integrators[1].getPotentialMaster(), space);
        sim.integrators[0].getMoveManager().addMCMove(wiggle[0]);//reference system
        wiggle[1] = new MCMoveClusterWiggleAceticAcid(sim,sim.integrators[1].getPotentialMaster(), space);
        sim.integrators[1].getMoveManager().addMCMove(wiggle[1]);//target system

        Box referenceBox = sim.box[0];
        Box targetBox = sim.box[1];
             
        if (false) {
            referenceBox.getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            targetBox.getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ColorSchemeByType colorScheme0 = (ColorSchemeByType) simGraphic.getDisplayBox(referenceBox).getColorScheme();
            ColorSchemeByType colorScheme1 = (ColorSchemeByType) simGraphic.getDisplayBox(targetBox).getColorScheme();
            DiameterHashByType diameterScheme0 = (DiameterHashByType) simGraphic.getDisplayBox(referenceBox).getDiameterHash();
            DiameterHashByType diameterScheme1 = (DiameterHashByType) simGraphic.getDisplayBox(targetBox).getDiameterHash();

            AtomType typeCH3 = species.getCH3Type();
            AtomType typeC = species.getCType();
            AtomType typeDBO = species.getDBOType();
            AtomType typeSBO = species.getSBOType();
            AtomType typeH = species.getHType();
            colorScheme0.setColor(typeCH3, Color.GREEN);
            diameterScheme0.setDiameter(typeCH3, 2 * 1.7);
            colorScheme0.setColor(typeC, Color.BLUE);
            colorScheme0.setColor(typeDBO, Color.RED);
            colorScheme0.setColor(typeSBO, Color.YELLOW);
            colorScheme0.setColor(typeH, Color.WHITE);
            colorScheme1.setColor(typeCH3, Color.GREEN);
            diameterScheme1.setDiameter(typeCH3, 2 * 1.7);
            colorScheme1.setColor(typeC, Color.BLUE);
            colorScheme1.setColor(typeDBO, Color.RED);
            colorScheme1.setColor(typeSBO, Color.YELLOW);
            colorScheme1.setColor(typeH, Color.WHITE);

            simGraphic.getDisplayBox(referenceBox).setShowBoundary(false);
            simGraphic.getDisplayBox(targetBox).setShowBoundary(false);
            simGraphic.makeAndDisplayFrame();
            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(targetBox).canvas).setBackgroundColor(Color.WHITE);
            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 100, false);
            sim.equilibrate(null, 200);
            sim.getController2().addActivity(new ActivityIntegrate2(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nBody+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize((int)steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));
        
        IAction progressReport = new IAction() {
            public void actionPerformed() {
            	System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nBody]+", error: "+ratioAndError[1]*HSB[nBody]);
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval(100);
        sim.integratorOS.getEventManager().addListener(progressReportListener);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(1000);
        sim.getController().actionPerformed();
        
        System.out.println("ideal reference step frequency "+sim.integratorOS.getIdealRefStepFraction());//optimize the uncertainty
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());//actually happened
      
        sim.printResults(HSB[nBody]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

    	public int nBody = 3;        
        public double temperature = 400;   // Kelvin        
        public long numSteps = 20000000;
        public double sigmaHSRef = 5.0;
    }

    
}


