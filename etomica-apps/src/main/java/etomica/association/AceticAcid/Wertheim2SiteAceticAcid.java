/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association.AceticAcid;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.association.BiasVolume2SiteAceticAcid;
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
import etomica.models.OPLS.P2AceticAcidTwoSite;
import etomica.models.OPLS.P2HardAssociationRefAceticAcidTwoSite;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.awt.*;

/**
 * Wertheim coefficients calculation for acetic acid using IMPROVED OPLS united-atom model.
 * For the calculation of diagrams having fR only, Hard-Sphere system is used as reference.
 * For the calculation of diagrams having F, Hard-Chain system is used as reference.  
 *   
 *@author Hye Min Kim
 * Fab, 2012
 *  
 */

public class Wertheim2SiteAceticAcid { 
	
    public static void main(String[] args) {
        WertheimParam params = new WertheimParam();
        Space space = Space3D.getInstance();
        ParseArgs.doParseArgs(params, args);

        final int nBody = params.nBody;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef=params.sigmaHSRef;
        double associationEnergy = params.associationEnergy;
        int numDiagram = params.numDiagram;
        boolean specialRef = params.specialRef;

        final double[] HSB = new double[9];

        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);

        System.out.println("IMPROVED OPLS Acetic acid with angle restriction Wertheim D"+numDiagram+" at "+temperature+"K");
        System.out.println("repulsive potential:>- "+associationEnergy+"cal/mole");
        double temperatureK = temperature;
        double associationEnergyCal = associationEnergy;
        temperature = Kelvin.UNIT.toSim(temperature);//T * kB
        Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
        associationEnergy = calPerMole.toSim(associationEnergy);

        //geometric association restriction (Consistent with NPT simulation)
        double innerOH = 1.8;
        double outerOH = 2.3;
        double maxCosOHO = -0.7;
        double maxCosHOC = -0.4;
        double minCosHOC = -0.9;

        System.out.println("innerOH "+innerOH+" outerOH "+outerOH+" maxCosOHO "+maxCosOHO+" maxCosHOC "+maxCosHOC+" minCosHOC "+minCosHOC);

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);

        PotentialGroup p = new PotentialGroup(2);

        P2AceticAcidTwoSite pR = new P2AceticAcidTwoSite(associationEnergy,false,3);//repulsive potential
        P2AceticAcidTwoSite pAB = new P2AceticAcidTwoSite(associationEnergy,true,0);//A-B association potential
        P2AceticAcidTwoSite pBA = new P2AceticAcidTwoSite(associationEnergy,true,1);//B-A association potential
        P2HardAssociationRefAceticAcidTwoSite pARef = new P2HardAssociationRefAceticAcidTwoSite();//hard-chain reference potential

        MayerEGeneral eR = new MayerEGeneral(pR);//repulsion eR function
	    MayerGeneral fR = new MayerGeneral(pR);//Mayer f function of reference part
	    MayerGeneral fAB = new MayerGeneral(pAB);
	    MayerGeneral fBA = new MayerGeneral(pBA);
	    MayerFunctionProductGeneral FAB = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fAB}, new double[]{1,1});
	    MayerFunctionProductGeneral FABFBA = new MayerFunctionProductGeneral(space, new MayerFunction[]{FAB,fBA}, new double[]{1,1});
	    MayerGeneral fARef = new MayerGeneral(pARef);

        int nBondTypes = 4;//fR, FAB, FAB2, eR
        ClusterBonds[] clusters = new ClusterBonds[0];
        int[][][] bondList = new int[nBondTypes][][];
        int [][][]refBondList = new int[nBondTypes][][];
        ClusterAbstract targetCluster = null;

        //for hard-chain reference value
        double hardShellVol = 4.0/3.0*Math.PI*(outerOH*outerOH*outerOH-innerOH*innerOH*innerOH);
        double hardSphereVol = 4.0/3.0*Math.PI*sigmaHSRef*sigmaHSRef*sigmaHSRef;
        double OHOFraction = 0.5*(1+maxCosOHO);
        double HOCFraction = 0.5*(maxCosHOC-minCosHOC);

        if (numDiagram >10 && numDiagram < 42 && !specialRef){
        	throw new RuntimeException("specialRef should be true.");
        }

        if (numDiagram == 2 || numDiagram == 3){//2 point diagrams having FAB or FABFBA
        	HSB[2] = -hardShellVol*OHOFraction*HOCFraction;
        	refBondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 5){//3 point diagrams having 1 FAB
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[3] = -hardSphereVol*253.2825;
        			else if  (temperatureK == 450) HSB[3] = -hardSphereVol*95.97541;
        			else if  (temperatureK == 500) HSB[3] = -hardSphereVol*46.83418;
        			else if  (temperatureK == 550) HSB[3] = -hardSphereVol*26.35662;
        			else if  (temperatureK == 600) HSB[3] = -hardSphereVol*16.87148;
        			else if  (temperatureK == 650) HSB[3] = -hardSphereVol*11.65923;
        		}
            	refBondList[0] = new int [][]{{1,2}};
            	refBondList[1] = new int [][]{{0,1}};
        	}
        	else {
        		HSB[3] = hardSphereVol*hardShellVol*OHOFraction*HOCFraction;
        		refBondList[0] = new int [][]{{1,2}};
        		refBondList[1] = new int [][]{{0,1}};
        	}
        }
        else if (numDiagram == 6){//3 point diagrams having 1 FABFBA
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[3] = -Math.pow(hardSphereVol,1)*173.8722;
        			else if  (temperatureK == 450) HSB[3] = -Math.pow(hardSphereVol,1)*48.31472;
        			else if  (temperatureK == 500) HSB[3] = -Math.pow(hardSphereVol,1)*17.67929;
        			else if  (temperatureK == 550) HSB[3] = -Math.pow(hardSphereVol,1)*7.784599;
        			else if  (temperatureK == 600) HSB[3] = -Math.pow(hardSphereVol,1)*3.819302;
        			else if  (temperatureK == 650) HSB[3] = -Math.pow(hardSphereVol,1)*2.100071;
        		}
        		refBondList[0] = new int [][]{{1,2}};
            	refBondList[2] = new int [][]{{0,1}};
        	}
        	else {
        		HSB[3] = hardSphereVol*hardShellVol*OHOFraction*HOCFraction;
        		refBondList[0] = new int [][]{{1,2}};
        		refBondList[1] = new int [][]{{0,1}};
        	}
        }
        else if (numDiagram == 7){//3 point diagrams having 2 FAB
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[3] = 253.2825*253.2825;
        			else if  (temperatureK == 450) HSB[3] = 95.97541*95.97541;
        			else if  (temperatureK == 500) HSB[3] = 46.83418*46.83418;
        			else if  (temperatureK == 550) HSB[3] = 26.35662*26.35662;
        			else if  (temperatureK == 600) HSB[3] = 16.87148*16.87148;
        			else if  (temperatureK == 650) HSB[3] = 11.65923*11.65923;
        		}
        		refBondList[1] = new int [][]{{0,1},{1,2}};
        	}
        	else {
        		HSB[3] = Math.pow(hardShellVol,2)*Math.pow(OHOFraction, 2)*Math.pow(HOCFraction,2);
        		refBondList[1] = new int [][]{{0,1},{1,2}};
        	}
        }
        else if (numDiagram >10 && numDiagram < 16){//4 point diagrams having 1 FAB
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = Math.pow(hardSphereVol,2)*253.2825;
        			else if  (temperatureK == 450) HSB[4] = Math.pow(hardSphereVol,2)*95.97541;
        			else if  (temperatureK == 500) HSB[4] = Math.pow(hardSphereVol,2)*46.83418;
        			else if  (temperatureK == 550) HSB[4] = Math.pow(hardSphereVol,2)*26.35662;
        			else if  (temperatureK == 600) HSB[4] = Math.pow(hardSphereVol,2)*16.87148;
        			else if  (temperatureK == 650) HSB[4] = Math.pow(hardSphereVol,2)*11.65923;
        		}
            	refBondList[0] = new int [][]{{1,2},{2,3}};
            	refBondList[1] = new int [][]{{0,1}};
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,2)*hardShellVol*OHOFraction*HOCFraction;
        		refBondList[0] = new int [][]{{1,2},{2,3}};
        		refBondList[1] = new int [][]{{0,1}};
        	}
        }
        else if (numDiagram >15 && numDiagram < 21){//4 point diagrams having 1 FABFBA
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = Math.pow(hardSphereVol,2)*173.8722;
        			else if  (temperatureK == 450) HSB[4] = Math.pow(hardSphereVol,2)*48.31472;
        			else if  (temperatureK == 500) HSB[4] = Math.pow(hardSphereVol,2)*17.67929;
        			else if  (temperatureK == 550) HSB[4] = Math.pow(hardSphereVol,2)*7.784599;
        			else if  (temperatureK == 600) HSB[4] = Math.pow(hardSphereVol,2)*3.819302;
        			else if  (temperatureK == 650) HSB[4] = Math.pow(hardSphereVol,2)*2.100071;
        		}
            	refBondList[0] = new int [][]{{1,2},{2,3}};
            	refBondList[2] = new int [][]{{0,1}};
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,2)*hardShellVol*OHOFraction*HOCFraction;
        		refBondList[0] = new int [][]{{1,2},{2,3}};
        		refBondList[1] = new int [][]{{0,1}};
        	}
        }
        else if (numDiagram >20 && numDiagram < 25){//4 point diagrams having 2 FAB
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = -Math.pow(hardSphereVol,1)*253.2825*253.2825;
        			else if  (temperatureK == 450) HSB[4] = -Math.pow(hardSphereVol,1)*95.97541*95.97541;
        			else if  (temperatureK == 500) HSB[4] = -Math.pow(hardSphereVol,1)*46.83418*46.83418;
        			else if  (temperatureK == 550) HSB[4] = -Math.pow(hardSphereVol,1)*26.35662*26.35662;
        			else if  (temperatureK == 600) HSB[4] = -Math.pow(hardSphereVol,1)*16.87148*16.87148;
        			else if  (temperatureK == 650) HSB[4] = -Math.pow(hardSphereVol,1)*11.65923*11.65923;
        		}
            	refBondList[0] = new int [][]{{1,2}};
            	refBondList[1] = new int [][]{{0,1},{2,3}};
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,1)*Math.pow((hardShellVol*OHOFraction*HOCFraction),2);
        		refBondList[0] = new int [][]{{1,2}};
        		refBondList[1] = new int [][]{{0,1},{2,3}};
        	}
        }
        else if (numDiagram >24 && numDiagram < 29){//4 point diagrams having 1 FAB and 1 FABFBA
        	if (specialRef){
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = -Math.pow(hardSphereVol,1)*253.2825*173.8722;
        			else if  (temperatureK == 450) HSB[4] = -Math.pow(hardSphereVol,1)*95.97541*48.31472;
        			else if  (temperatureK == 500) HSB[4] = -Math.pow(hardSphereVol,1)*46.83418*17.67929;
        			else if  (temperatureK == 550) HSB[4] = -Math.pow(hardSphereVol,1)*26.35662*7.784599;
        			else if  (temperatureK == 600) HSB[4] = -Math.pow(hardSphereVol,1)*16.87148*3.819302;
        			else if  (temperatureK == 650) HSB[4] = -Math.pow(hardSphereVol,1)*11.65923*2.100071;
        		}
            	refBondList[0] = new int [][]{{1,2}};
            	refBondList[1] = new int [][]{{0,1}};
            	refBondList[2] = new int [][]{{2,3}};
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,1)*Math.pow((hardShellVol*OHOFraction*HOCFraction),2);
            	refBondList[0] = new int [][]{{1,2}};
            	refBondList[1] = new int [][]{{0,1},{2,3}};
        	}
        }
        else if ((numDiagram >28 && numDiagram < 32)){//4 point diagrams having 2 FABFBA
        	if (specialRef){//D3*D3*sigmaHSRef is the reference value
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = -Math.pow(hardSphereVol,1)*173.8722*173.8722;
        			else if  (temperatureK == 450) HSB[4] = -Math.pow(hardSphereVol,1)*48.31472*48.31472;
        			else if  (temperatureK == 500) HSB[4] = -Math.pow(hardSphereVol,1)*17.67929*17.67929;
        			else if  (temperatureK == 550) HSB[4] = -Math.pow(hardSphereVol,1)*7.784599*7.784599;
        			else if  (temperatureK == 600) HSB[4] = -Math.pow(hardSphereVol,1)*3.819302*3.819302;
        			else if  (temperatureK == 650) HSB[4] = -Math.pow(hardSphereVol,1)*2.100071*2.100071;
        		}
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,1)*Math.pow((hardShellVol*OHOFraction*HOCFraction),2);
        	}
        	refBondList[0] = new int [][]{{1,2}};
        	refBondList[1] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram >31 && numDiagram < 37){//4 point diagrams having 2 adjacent FAB
        	if (specialRef){//D2*D2*sigmaHSRef is the reference value
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = -Math.pow(hardSphereVol,1)*253.2825*253.2825;
        			else if  (temperatureK == 450) HSB[4] = -Math.pow(hardSphereVol,1)*95.97541*95.97541;
        			else if  (temperatureK == 500) HSB[4] = -Math.pow(hardSphereVol,1)*46.83418*46.83418;
        			else if  (temperatureK == 550) HSB[4] = -Math.pow(hardSphereVol,1)*26.35662*26.35662;
        			else if  (temperatureK == 600) HSB[4] = -Math.pow(hardSphereVol,1)*16.87148*16.87148;
        			else if  (temperatureK == 650) HSB[4] = -Math.pow(hardSphereVol,1)*11.65923*11.65923;
        		}
        	}
        	else {
        		HSB[4] = -Math.pow(hardSphereVol,1)*Math.pow((hardShellVol*OHOFraction*HOCFraction),2);
        	}
        	refBondList[0] = new int [][]{{2,3}};
        	refBondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram >36 && numDiagram < 42){//4 point diagrams having 3 adjacent FAB
        	if (specialRef){//D2*D2*D2 is the reference value
        		if (associationEnergyCal == 5000.0){
        			if (temperatureK == 400) HSB[4] = 253.2825*253.2825*253.2825;
        			else if  (temperatureK == 450) HSB[4] = 95.97541*95.97541*95.97541;
        			else if  (temperatureK == 500) HSB[4] = 46.83418*46.83418*46.83418;
        			else if  (temperatureK == 550) HSB[4] = 26.35662*26.35662*26.35662;
        			else if  (temperatureK == 600) HSB[4] = 16.87148*16.87148*16.87148;
        			else if  (temperatureK == 650) HSB[4] = 11.65923*11.65923*11.65923;
        		}
        	}
        	else {
        		HSB[4] = -Math.pow((hardShellVol*OHOFraction*HOCFraction),3);
        	}
        	refBondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        if (numDiagram == 1) {
        	bondList[0]=new int [][]{{0,1}};
        }
        else if (numDiagram == 2) {
        	bondList[1]=new int [][]{{0,1}};
        }
        else if (numDiagram == 3) {
        	bondList[2]=new int [][]{{0,1}};
        }
        else if (numDiagram == 4) {
        	bondList[0]=new int [][]{{0,1},{1,2},{2,0}};
        }
        else if (numDiagram == 5) {
        	bondList[0] = new int [][]{{1,2},{2,0}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 6) {
        	bondList[0] = new int [][]{{1,2},{2,0}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 7) {
        	bondList[0] = new int [][]{{2,0}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 8) {
        	bondList[0]=new int [][]{{0,1},{1,2},{2,3},{3,0}};
        }
        else if (numDiagram == 9) {
        	bondList[0]=new int [][]{{0,1},{1,2},{2,3},{3,0},{0,2}};
        }
        else if (numDiagram == 10) {
        	bondList[0]=new int [][]{{0,1},{1,2},{2,3},{3,0},{0,2},{3,1}};
        }
        else if (numDiagram == 11) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 12) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{0,2}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 13) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{1,3}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 14) {
        	bondList[0]=new int [][]{{1,2},{3,0},{0,2},{1,3}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 15) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{0,2},{1,3}};
        	bondList[1] = new int [][]{{0,1}};
        }
        else if (numDiagram == 16) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 17) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{0,2}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 18) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{1,3}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 19) {
        	bondList[0]=new int [][]{{1,2},{3,0},{0,2},{1,3}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 20) {
        	bondList[0]=new int [][]{{1,2},{2,3},{3,0},{0,2},{1,3}};
        	bondList[2] = new int [][]{{0,1}};
        }
        else if (numDiagram == 21) {
        	bondList[0]=new int [][]{{1,2},{3,0}};
        	bondList[1] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 22) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0}};
        	bondList[1] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 23) {
        	bondList[0]=new int [][]{{1,2},{3,0},{1,3}};
        	bondList[1] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 24) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0},{1,3}};
        	bondList[1] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 25) {
        	bondList[0]=new int [][]{{1,2},{3,0}};
        	bondList[1] = new int [][]{{0,1}};
        	bondList[2] = new int [][]{{2,3}};
        }
        else if (numDiagram == 26) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0}};
        	bondList[1] = new int [][]{{0,1}};
        	bondList[2] = new int [][]{{2,3}};
        }
        else if (numDiagram == 27) {
        	bondList[0]=new int [][]{{1,2},{3,0},{1,3}};
        	bondList[1] = new int [][]{{0,1}};
        	bondList[2] = new int [][]{{2,3}};
        }
        else if (numDiagram == 28) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0},{1,3}};
        	bondList[1] = new int [][]{{0,1}};
        	bondList[2] = new int [][]{{2,3}};
        }
        else if (numDiagram == 29) {
        	bondList[0]=new int [][]{{1,2},{3,0}};
        	bondList[2] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 30) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0}};
        	bondList[2] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 31) {
        	bondList[0]=new int [][]{{1,2},{3,0},{2,0},{1,3}};
        	bondList[2] = new int [][]{{0,1},{2,3}};
        }
        else if (numDiagram == 32) {
        	bondList[0]=new int [][]{{2,3},{3,0}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 33) {
        	bondList[0]=new int [][]{{2,3},{3,0},{0,2}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 34) {
        	bondList[0]=new int [][]{{2,3},{3,1},{0,2}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 35) {
        	bondList[0]=new int [][]{{2,3},{3,0},{1,3}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 36) {
        	bondList[0]=new int [][]{{2,3},{3,0},{0,2},{1,3}};
        	bondList[1] = new int [][]{{0,1},{1,2}};
        }
        else if (numDiagram == 37) {
        	bondList[0]=new int [][]{{3,0}};
        	bondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        else if (numDiagram == 38) {
        	bondList[0]=new int [][]{{3,0},{0,2}};
        	bondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        else if (numDiagram == 39) {
        	bondList[0]=new int [][]{{2,0},{1,3}};
        	bondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        else if (numDiagram == 40) {
        	bondList[0]=new int [][]{{3,0},{1,3}};
        	bondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        else if (numDiagram == 41) {
        	bondList[0]=new int [][]{{3,0},{0,2},{1,3}};
        	bondList[1] = new int [][]{{0,1},{1,2},{2,3}};
        }
        else {
        	throw new RuntimeException("This is strange");
        }

        clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
        targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FAB,FABFBA,eR, fARef});
        if (numDiagram == 1) {
        	System.out.println("Flipping is applied");
        	((ClusterSum)targetCluster).setCaching(false);
        	targetCluster = new ClusterCoupledFlipped(targetCluster, space);
        }
        targetCluster.setTemperature(temperature);

        ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody>3, eRef, true);
        if (numDiagram == 2 || numDiagram == 3 || (numDiagram > 4 && numDiagram < 8)||(numDiagram > 10 && numDiagram < 25)) {
    		ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
        	if (specialRef){
        		System.out.println("D2 (or D3) * Hard-Sphere reference is used.");
            	refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, FAB, FABFBA});
        	}
        	else {
        		System.out.println("Hard Chain Reference is used.");

        		refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, fARef});
        	}
        }
        else if ((numDiagram > 24 && numDiagram < 29)) {
        	ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
        	System.out.println("D2 * D3 * Hard-Sphere reference is used.");
        	refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, FAB, FABFBA});
        }
        else if ((numDiagram > 28 && numDiagram < 32)) {
        	ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
        	System.out.println("D3 * D3 * Hard-Sphere reference is used.");
        	refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, FABFBA});
        }
        else if ((numDiagram > 31 && numDiagram < 37)) {
        	ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
        	System.out.println("D2 * D2 * Hard-Sphere reference is used.");
        	refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, FAB});
        }
        else if ((numDiagram > 36 && numDiagram < 42)) {
        	ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
        	System.out.println("D2 * D2 * D2 is used.");
        	refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fRef, FAB});
        }
        refCluster.setTemperature(temperature);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2 (space,new SpeciesAceticAcid(space),temperature,refCluster,targetCluster, false);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        SpeciesAceticAcid species = (SpeciesAceticAcid)sim.getSpecies(0);

        AceticAcidModPotentialHelper.initPotential(space, species, p);
        AceticAcidModPotentialHelper.initPotential(space, species, pR);
        AceticAcidModPotentialHelper.initPotential(space, species, pAB);
        AceticAcidModPotentialHelper.initPotential(space, species, pBA);

        BiasVolume2SiteAceticAcid bv = new BiasVolume2SiteAceticAcid(space, sim.getRandom());

        pR.setBiasVolume(bv);
        pAB.setBiasVolume(bv);
        pBA.setBiasVolume(bv);
        pARef.setBiasVolume(bv);

        System.out.println("B"+nBody+"HS: "+HSB[nBody]);

       //set the initial configuration to satisfy the associations among particles
        if (numDiagram == 2 || numDiagram == 5 || (numDiagram >10 && numDiagram < 16)){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FAB);
        	configuration.initializeCoordinates(sim.box[1],false,false);
        	configuration.initializeCoordinates(sim.box[0],false,false);
        }
        else if (numDiagram == 3 || numDiagram == 6|| (numDiagram >15 && numDiagram < 21)){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FABFBA);
        	configuration.initializeCoordinates(sim.box[1],false,false);
        	configuration.initializeCoordinates(sim.box[0],false,false);
        }
        else if (numDiagram == 7 || (numDiagram >31 && numDiagram < 37)){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FAB);
        	configuration.initializeCoordinates(sim.box[1],true,false);
        	configuration.initializeCoordinates(sim.box[0],true,false);
        }
        else if (numDiagram >20 && numDiagram < 25){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FAB);
        	configuration.initializeCoordinates(sim.box[1],false,true);
        	configuration.initializeCoordinates(sim.box[0],false,true);
        }
        else if (numDiagram >24 && numDiagram < 29){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FAB, FABFBA);
        	configuration.initializeCoordinates2(sim.box[1],false,true);
        	configuration.initializeCoordinates2(sim.box[0],false,true);
        }
        else if (numDiagram >28 && numDiagram < 32){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FABFBA);
        	configuration.initializeCoordinates(sim.box[1],false,true);
        	configuration.initializeCoordinates(sim.box[0],false,true);
        }
        else if (numDiagram >36 && numDiagram < 42){
        	ConfigurationClusterAceticAcid configuration = new ConfigurationClusterAceticAcid(space, sim.getRandom(), FAB);
        	configuration.initializeCoordinates(sim.box[1],true,true);
        	configuration.initializeCoordinates(sim.box[0],true,true);
        }

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

        //utorsional = c1*[1+cos(phi+f1)]+c2*[1-cos(2*phi)], c1/kB = 630K, c2/kB = 1562.4K, f1=180째(for OCOH), 0째(for CCOH)

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

        if (true) {
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
            diameterScheme0.setDiameter(typeCH3, 2*1.7);
            colorScheme0.setColor(typeC, Color.BLUE);
            colorScheme0.setColor(typeDBO, Color.RED);
            colorScheme0.setColor(typeSBO, Color.YELLOW);
            colorScheme0.setColor(typeH, Color.WHITE);
            colorScheme1.setColor(typeCH3, Color.GREEN);
            diameterScheme1.setDiameter(typeCH3, 2*1.7);
            colorScheme1.setColor(typeC, Color.BLUE);
            colorScheme1.setColor(typeDBO, Color.RED);
            colorScheme1.setColor(typeSBO, Color.YELLOW);
            colorScheme1.setColor(typeH, Color.WHITE);

            simGraphic.getDisplayBox(referenceBox).setShowBoundary(false);
            simGraphic.getDisplayBox(targetBox).setShowBoundary(false);
            simGraphic.makeAndDisplayFrame();
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(targetBox).canvas).setBackgroundColor(Color.WHITE);
            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 100, false);
            sim.equilibrate(null, 200);
            sim.getController().addActivity(new ActivityIntegrate2(sim.integratorOS));

            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+numDiagram+"_"+temperature : null;
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
        		System.out.println("abs average: "+ratioAndError[0]*HSB[nBody]+", error: "+ratioAndError[1]*Math.abs(HSB[nBody]));
        	}
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval(100);
        sim.integratorOS.getEventManager().addListener(progressReportListener);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), 1000);

        System.out.println("ideal reference step frequency "+sim.integratorOS.getIdealRefStepFraction());//optimize the uncertainty
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());//actually happened

        sim.printResults(HSB[nBody]);
 }

    /**
     * Inner class for parameters
     */
    public static class WertheimParam extends ParameterBase {

    	public int nBody = 3;       
    	public double temperature = 400.0;// Kelvin       
    	public long numSteps = 10000000;
    	public double sigmaHSRef = 5.0;
    	public double associationEnergy = 5000.0;
    	public int numDiagram = 4;
    	public boolean specialRef = false;
    }
    	
   
}


