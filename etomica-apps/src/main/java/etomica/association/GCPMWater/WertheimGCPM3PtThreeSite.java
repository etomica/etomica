/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association.GCPMWater;

import etomica.action.IAction;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.listener.IntegratorListenerAction;
import etomica.models.water.P2HardAssociationGCPMReference;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.models.water.SpeciesWater4P;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap;

import java.awt.*;

/**
 * repulsive potential: energy of pair is greater than -association energy(cal/mol)
 * attractive potential: energy of pair is less than -association energy(cal/mol)
 * Wertheim's three attraction-site model of GCPM 
 * Three point diagram
 * Three Site model
 *
 * @author Hye Min Kim
 */

public class WertheimGCPM3PtThreeSite {

	public static void main(String[] args) {
//    	Unit calPerMoles = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
//    	System.out.println(calPerMoles.fromSim(Kelvin.UNIT.toSim(110)));
//    	System.exit(1);
		VirialAssociatingFluidParam params = new VirialAssociatingFluidParam();
		Space space = Space3D.getInstance();
		
		int numDiagram=params.numDiagram;
		int diagramIndex = params.diagramIndex;
		final int nBody = 3;//W3
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
        System.out.println("Wertheim's three association site model");
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
		
		int nBondTypes = 5;//fR,FAC,FCA,eR
		ClusterBonds[] clusters = new ClusterBonds[0];
		int[][][] bondList = new int[nBondTypes][][];
		int [][][]refBondList = new int[nBondTypes][][];
		ClusterSum targetCluster = null;

		System.out.println("Digaram "+numDiagram+"-"+diagramIndex);
		if (numDiagram == 3 ) {
			bondList[0]=new int [][]{{0,1},{1,2},{2,0}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FAC,eR});
			
		}  else if (numDiagram == 6){
			if (diagramIndex == 1) {
				HSB[3] = -35.238*-35.238;//This value is from direct sampling
				refBondList[1] = new int [][]{{0,1}};
				refBondList[0] = new int [][]{{1,2}};
				bondList[2] = new int [][]{{0,1}};
				bondList[1] = new int [][]{{1,2}};
				bondList[3] = new int [][]{{2,0}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,eR});
			}
			else if (diagramIndex == 2) {
				HSB[3] = -35.238*-35.238;//This value is from direct sampling
				refBondList[0] = new int [][]{{0,1}};
				refBondList[1] = new int [][]{{1,2}};
				bondList[1] = new int [][]{{0,1}};
				bondList[2] = new int [][]{{1,2}};
				bondList[3] = new int [][]{{2,0}};
				clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
				targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,eR});
			}
		}  	else if (numDiagram == 7){
			if (diagramIndex == 1){

			HSB[3] = -35.238*-35.238;//This value is from direct sampling
			refBondList[0] = new int [][]{{0,1},{1,2}};
			bondList[1] = new int [][]{{0,1},{1,2}};
			bondList[2] = new int [][]{{2,0}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,eR});
		}
		else if (diagramIndex == 2){
			HSB[3] = -35.238*-35.238;//This value is from direct sampling
			refBondList[0] = new int [][]{{0,1}};
			refBondList[1] = new int [][]{{1,2}};
			bondList[1] = new int [][]{{0,1},{0,2}};
			bondList[2] = new int [][]{{1,2}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,eR});
		}
		}
		else if (numDiagram ==4){
			bondList[0] = new int [][]{{1,2},{2,0}};
			bondList[1] = new int [][]{{0,1}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FAC,eR});
		}	else if (numDiagram ==5){
			if (diagramIndex == 1) {
				HSB[3] = -35.238*-35.238;//This value is from direct sampling
				refBondList[0] = new int [][]{{0,1},{1,2}};
				bondList[0] = new int [][]{{2,0}};
				bondList[1] = new int [][]{{0,1},{1,2}};
				clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
				targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,eR});
			}	else if (diagramIndex == 2) {
				HSB[3] = -35.238*-35.238;//This value is from direct sampling
				refBondList[0] = new int [][]{{0,1}};
				refBondList[2] = new int [][]{{1,2}};
				bondList[0] = new int [][]{{2,0}};
				bondList[1] = new int [][]{{0,1}};
				bondList[3] = new int [][]{{1,2}};
				clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
				targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,FBC,eR});
			}
		}	
			else {
			throw new RuntimeException("This is strange");
		}

		System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
		ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody>3, eRef, true);
		if ((numDiagram == 5 ||numDiagram == 6 ||numDiagram == 7)){
			ClusterBonds refBonds = new ClusterBonds(4,refBondList, false);
			refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fCARef4mer, fACRef4mer, fBCRef4mer});
		}
        refCluster.setTemperature(temperature);
        targetCluster.setTemperature(temperature);
		final SimulationVirialOverlap sim = new SimulationVirialOverlap(space, new SpeciesFactoryWaterGCPM(), temperature,refCluster,targetCluster);
		//ConfigurationClusterMove configuration = new ConfigurationClusterMove(space, sim.getRandom());
		ConfigurationClusterWertheimGCPM configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pAC);
		if (numDiagram == 3 || numDiagram == 4 ) {
			configuration.initializeCoordinates(sim.box[1]);
		}
		if (numDiagram ==5 &&diagramIndex==1){//diagram 5-1
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
			if (numDiagram ==7 &&diagramIndex==1){//diagram 7-1
				configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pCA, pCA, pCA);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates5(sim.box[1]);	
			}
			if (numDiagram ==7 &&diagramIndex==2){//diagram 7-2
				configuration = new ConfigurationClusterWertheimGCPM(space, sim.getRandom(), pCA, pAC, pCA);
				configuration.initializeCoordinates3(sim.box[0]);	
				configuration.initializeCoordinates5(sim.box[1]);		
		}	
		sim.setAccumulatorBlockSize((int)numSteps*10);
		sim.integratorOS.setNumSubSteps(1000);	

		if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            SpeciesWater4P species = (SpeciesWater4P)sim.getSpecies(0);
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
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null,20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
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
                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nBody]+", error: "+ratioAndError[1]*HSB[nBody]);
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval((int)(numSteps/10));
        sim.integratorOS.getEventManager().addListener(progressReportListener);
        
        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());

		sim.printResults(HSB[nBody]);
	}


	public static class VirialAssociatingFluidParam extends ParameterBase {
		public double temperature = 630;//reduced temperature
		public double sigmaHSRef = 3.2;
		public long numSteps = 5000;
		public double associationEnergy = 2500.0;
		public int numDiagram = 6;
		public int diagramIndex = 2;
		
	}

}
