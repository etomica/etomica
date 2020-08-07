/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association.GCPMWater;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.P2HardAssociationGCPMReference;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.models.water.SpeciesWater4P;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.awt.*;

/**
 * repulsive potential: energy of pair is greater than -3000cal/mol
 * attractive potential: energy of pair is less than -3000cal/mol
 * Wertheim's single attraction-site model of GCPM 
 * Three Site model
 * e-bond in 3 body term is decomposed into eR and F
 *
 * @author Hye Min Kim
 */

public class WertheimGCPM3PtEBondDecomp {

	public static void main(String[] args) {
//    	Unit calPerMoles = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
//    	System.out.println(calPerMoles.fromSim(Kelvin.UNIT.toSim(110)));
//    	System.exit(1);
		VirialAssociatingFluidParam params = new VirialAssociatingFluidParam();
		Space space = Space3D.getInstance();
		
		int numDiagram=params.numDiagram;
		int diagramIndex=params.diagramIndex;
		final int nBody = 3;//B3
		double temperature=params.temperature;
		double sigmaHSRef = params.sigmaHSRef;
		double associationEnergy = params.associationEnergy;
		long numSteps = params.numSteps;
		
		if (args.length == 6) {
        	temperature = Double.parseDouble(args[0]);
            sigmaHSRef = Double.parseDouble(args[1]);
            numSteps = Integer.parseInt(args[2]);
            associationEnergy = Double.parseDouble(args[3]);
            numDiagram = Integer.parseInt(args[4]);
            diagramIndex = Integer.parseInt(args[5]);
            
        } else if (args.length != 0){
        	throw new IllegalArgumentException("Wrong number of arguments");
        }

		final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        System.out.println("Wertheim's three association site model && e-bond in 3 body term is decomposed into eR and F");
        System.out.println("repulsive potential: GCPM potential - negative charge electrostatic part, attractive potential: negative charge electrostatic potential");
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("temperature: "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println("association Energy: "+associationEnergy+"cal/mole");
    	Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	associationEnergy = calPerMole.toSim(associationEnergy);

		PNWaterGCPMThreeSite pR = new PNWaterGCPMThreeSite(space,associationEnergy,false); //GCPM potential
		PNWaterGCPMThreeSite pAC = new PNWaterGCPMThreeSite(space,associationEnergy,true); //GCPM potential
		PNWaterGCPMThreeSite pCA = new PNWaterGCPMThreeSite(space,associationEnergy,true); //GCPM potential
		PNWaterGCPMThreeSite pBC = new PNWaterGCPMThreeSite(space,associationEnergy,true); //GCPM potential
		pAC.setBondType(1);
		pCA.setBondType(3);
		pBC.setBondType(2);
		MayerEGeneral eR = new MayerEGeneral(pR);//repulsion eR function
		MayerGeneral fR = new MayerGeneral(pR);//Mayer f function of reference part
		MayerGeneral fAC = new MayerGeneral(pAC);//Mayer f function of association part
		MayerGeneral fCA = new MayerGeneral(pCA);//Mayer f function of association part
		MayerGeneral fBC = new MayerGeneral(pBC);//Mayer f function of association part
		MayerFunctionProductGeneral FAC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fAC}, new double[]{1,1});//F=eR*fA
		MayerFunctionProductGeneral FCA = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fCA}, new double[]{1,1});//F=eR*fA
		MayerFunctionProductGeneral FBC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fBC}, new double[]{1,1});//F=eR*fA
		
		MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
		P2HardAssociationGCPMReference pCARef = new P2HardAssociationGCPMReference(space, true);
		pCARef.setBondType(3);
		MayerGeneral fCARef4mer = new MayerGeneral(pCARef);
		P2HardAssociationGCPMReference pACRef = new P2HardAssociationGCPMReference(space, true);
		pACRef.setBondType(1);
		MayerGeneral fACRef4mer = new MayerGeneral(pACRef);
		P2HardAssociationGCPMReference pBCRef = new P2HardAssociationGCPMReference(space, true);
		pBCRef.setBondType(2);
		MayerGeneral fBCRef4mer = new MayerGeneral(pBCRef);
		
		int nBondTypes = 5;//fR,FAC,FCA,FBC,eR
		ClusterBonds[] clusters = new ClusterBonds[0];
		int[][][] bondList = new int[nBondTypes][][];
		int [][][]refBondList = new int[nBondTypes][][];
		ClusterSumPolarizableWertheimProduct targetCluster = null;

		System.out.println("Digaram "+numDiagram+"-"+diagramIndex+"E");
		
		if (numDiagram == 5 ||numDiagram == 6||numDiagram == 7){
			HSB[3] = -35.238*-35.238;//This value is from direct sampling
			System.out.println("use hard-chain reference && value = -35.238*-35.238");
		}
		
		if (numDiagram == 3) {
			bondList[4]=new int [][]{{0,1},{1,2},{2,0}};
			
		}  else if (numDiagram == 6){
			if (diagramIndex == 1) {
				refBondList[1] = new int [][]{{0,1}};
				refBondList[0] = new int [][]{{1,2}};
				bondList[1] = new int [][]{{0,1}};
				bondList[3] = new int [][]{{1,2}};
				bondList[4] = new int [][]{{2,0}};
			}
			else if (diagramIndex == 2) {
				refBondList[0] = new int [][]{{0,1}};
				refBondList[1] = new int [][]{{1,2}};
				bondList[3] = new int [][]{{0,1}};
				bondList[1] = new int [][]{{1,2}};
				bondList[4] = new int [][]{{2,0}};
			}
		}  	
		else if (numDiagram ==4){
			bondList[4] = new int [][]{{1,2},{2,0}};
			bondList[1] = new int [][]{{0,1}};
		}	else if (numDiagram ==5){
			if (diagramIndex == 1) {
				refBondList[0] = new int [][]{{0,1},{1,2}};
				bondList[4] = new int [][]{{2,0}};
				bondList[3] = new int [][]{{0,1},{1,2}};
			}	else if (diagramIndex == 2) {
				refBondList[0] = new int [][]{{0,1}};
				refBondList[2] = new int [][]{{1,2}};
				bondList[4] = new int [][]{{2,0}};
				bondList[3] = new int [][]{{0,1}};
				bondList[2] = new int [][]{{1,2}};
			}
		}
			else {
			throw new RuntimeException("This is strange");
		}

		System.out.println("B3HS: " + HSB[3] + " = " + (HSB[3] / (HSB[2] * HSB[2])) + " B2HS^2");
		ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody > 3, eRef, true);
		if ((numDiagram == 5 || numDiagram == 6 || numDiagram == 7)) {
			ClusterBonds refBonds = new ClusterBonds(4, refBondList, false);
			refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double[]{1}, new MayerFunction[]{fCARef4mer, fACRef4mer, fBCRef4mer});
		}
		refCluster.setTemperature(temperature);
		clusters = (ClusterBonds[]) Arrays.addObject(clusters, new ClusterBonds(nBody, bondList, false));
		targetCluster = new ClusterSumPolarizableWertheimProduct(clusters, new double[]{1}, new MayerFunction[]{fR, FAC, FBC, FCA, eR});
		targetCluster.setTemperature(temperature);
		SpeciesWater4P species = new SpeciesWater4P(space);
		species.setConformation(new ConformationWaterGCPM(space));
		final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, temperature, refCluster, targetCluster);
		ConfigurationClusterWertheimGCPM configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pAC);
		if (numDiagram == 3) {
			configuration.initializeCoordinatesER(sim.box[1]);
		}
		if (numDiagram == 4) {
			configuration.initializeCoordinates(sim.box[1]);
		}

		if (numDiagram == 5 && diagramIndex == 1) {//diagram 5-1
			configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pCA, pCA);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates3(sim.box[1]);	
			}
			if (numDiagram ==5 &&diagramIndex==2){//diagram 5-2
				configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pCA, pBC);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates3(sim.box[1]);	
			}
			if (numDiagram ==6 &&diagramIndex==2){//diagram 6-2
				configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pCA, pAC);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates3(sim.box[1]);	
			}
			if (numDiagram ==6 &&diagramIndex==1){//diagram 6-1
				configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pAC, pCA);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates3(sim.box[1]);	
			}
						
		sim.setAccumulatorBlockSize((int)numSteps*10);
		sim.integratorOS.setNumSubSteps(1000);	

		if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);

            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            simGraphic.makeAndDisplayFrame();
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box[1]).canvas).setBackgroundColor(Color.WHITE);
            sim.integratorOS.setNumSubSteps(100000);
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
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+numDiagram+"_"+diagramIndex+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        //sim.setRefPref(50.0);
        sim.initRefPref(refFileName, numSteps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, numSteps/20);    
        sim.setAccumulatorBlockSize((int)numSteps);
                
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize());

        IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
				double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nBody]+", error: "+ratioAndError[1]*HSB[nBody]);
			}
		};
		IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
		progressReportListener.setInterval((int) (numSteps / 10));
		sim.integratorOS.getEventManager().addListener(progressReportListener);

		sim.integratorOS.getMoveManager().setEquilibrating(false);
		sim.ai.setMaxSteps(numSteps);
		sim.getController().actionPerformed();

		System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
		System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

		sim.printResults(HSB[nBody]);
	}
    
	
	public static class VirialAssociatingFluidParam extends ParameterBase {
		public double temperature = 600;//reduced temperature
		public double sigmaHSRef = 5.0;
		public long numSteps = 10000;
		public double associationEnergy = 3000.0;
		public int numDiagram = 3;
		public int diagramIndex = 1;
		
	}

}
