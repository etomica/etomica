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
import etomica.models.water.*;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.awt.*;

/**
 * repulsive potential: energy of pair is greater than -2500cal/mol
 * attractive potential: energy of pair is less than -2500cal/mol
 * Wertheim's single attraction-site model of GCPM 
 * Three Site model
 * e-bond in 4 body term is decomposed into eR and F
 * final revision Aug 2011
 *
 * @author Hye Min Kim
 */

public class WertheimGCPM4PtEBondDecomp {

	public static void main(String[] args) {
		WertheimAssociatingFluidParam params = new WertheimAssociatingFluidParam();
		Space space = Space3D.getInstance();
    	ParseArgs.doParseArgs(params, args);
        double deltaCut = 100;
		
		int numDiagram=params.numDiagram;
		int diagramIndex=params.diagramIndex;
		final int nBody = 4;//W4
		double temperature=params.temperature;
		double sigmaHSRef = params.sigmaHSRef;
		double associationEnergy = params.associationEnergy;
		long numSteps = params.numSteps;
		boolean flip = params.flip;
		boolean bondingAngleRestriction = params.bondingAngleRestriction;

		final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        System.out.println("Non-additive Wertheim diagram, final version, fij in delta B4 is NOT decomposed into fR and fAssociation");
        System.out.println("Wertheim's three association site model && e-bond in 4 body term is decomposed into eR and F");
        System.out.println("sigmaHSRef: "+sigmaHSRef);
//        System.out.println("B2HS: "+HSB[2]);
//        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("temperature: "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println("association Energy: "+associationEnergy+"cal/mol");//association energy Unit: cal/mol
        if (bondingAngleRestriction){
        	System.out.println("This calculation includes bonding angle criteria(cosTheta1<-0.5,cosTheta3<0)");
        }
    	Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	associationEnergy = calPerMole.toSim(associationEnergy);
    	
        IPotentialMolecular pTarget = new PNWaterGCPM(space);
        
        MayerEGeneral e = new MayerEGeneral(pTarget);
        
		PNWaterGCPMThreeSite pR = new PNWaterGCPMThreeSite(space,associationEnergy,false, bondingAngleRestriction); //GCPM potential
		PNWaterGCPMThreeSite pAC = new PNWaterGCPMThreeSite(space,associationEnergy,true, bondingAngleRestriction); //assciation between site A on mol1 and C on mol2
		PNWaterGCPMThreeSite pCA = new PNWaterGCPMThreeSite(space,associationEnergy,true, bondingAngleRestriction); //assciation between site C on mol1 and A on mol2
		PNWaterGCPMThreeSite pBC = new PNWaterGCPMThreeSite(space,associationEnergy,true, bondingAngleRestriction); //assciation between site B on mol1 and C on mol2
		PNWaterGCPMThreeSite pCB = new PNWaterGCPMThreeSite(space,associationEnergy,true, bondingAngleRestriction); //assciation between site C on mol1 and B on mol2
		pAC.setBondType(1);
		pCA.setBondType(3);
		pBC.setBondType(2);
		pCB.setBondType(4);
		MayerEGeneral eR = new MayerEGeneral(pR);//repulsion eR function
		MayerGeneral fR = new MayerGeneral(pR);//Mayer f function of reference part
		MayerGeneral fAC = new MayerGeneral(pAC);//Mayer f function of association part
		MayerGeneral fCA = new MayerGeneral(pCA);//Mayer f function of association part
		MayerGeneral fBC = new MayerGeneral(pBC);//Mayer f function of association part
		MayerGeneral fCB = new MayerGeneral(pCB);//Mayer f function of association part
		MayerFunctionProductGeneral FAC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fAC}, new double[]{1,1});//association bond between site A on mol1 and C on mol2
		MayerFunctionProductGeneral FCA = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fCA}, new double[]{1,1});//association bond between site C on mol1 and A on mol2
		MayerFunctionProductGeneral FBC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fBC}, new double[]{1,1});//association bond between site B on mol1 and C on mol2
		MayerFunctionProductGeneral FCB = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fCB}, new double[]{1,1});//association bond between site C on mol1 and B on mol2
		
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
		
		int nBondTypes = 7;//fR,FCA,FAC,FCB,FBC,eR
		ClusterBonds[] clusters = new ClusterBonds[0];
		int[][][] bondList = new int[nBondTypes][][];
		int [][][]refBondList = new int[nBondTypes][][];
		ClusterAbstract targetCluster = null;
		
		System.out.println("Diagram "+numDiagram+"-"+diagramIndex+"E "+"flip "+flip);
		
		if (numDiagram == 12 ||numDiagram == 19||numDiagram == 26||numDiagram == 34){
			HSB[4] = -35.238*-35.238*-35.238;//This value is from direct sampling
			System.out.println("use hard-chain reference && value = -35.238*-35.238*-35.258 = 43755.611101272");
		}
            		
		if (numDiagram == 12){			
				if (diagramIndex == 5){
					refBondList[0] = new int [][]{{0,1},{2,3}};
					refBondList[1] = new int [][]{{1,2}};
					bondList[5]=new int [][]{{3,0},{0,2},{1,3}};
					bondList[1]=new int [][]{{0,1},{2,3}};
					bondList[2]=new int [][]{{1,2}};
				}
		}	else if (numDiagram == 16){
				if (diagramIndex == 3){
					bondList[5]=new int [][]{{2,3},{3,0},{1,3},{0,2}};
					bondList[2]=new int [][]{{0,1}};
					bondList[1]=new int [][]{{1,2}};
				}	else if (diagramIndex == 4){
					bondList[5]=new int [][]{{2,3},{3,0},{1,3},{0,2}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2}};
				}
		}	else if (numDiagram == 19){
				if (diagramIndex == 3){
					refBondList[1] = new int [][]{{0,1}};
					refBondList[2] = new int [][]{{2,3}};
					refBondList[0] = new int [][]{{1,2}};		
					bondList[5]=new int [][]{{3,0},{1,3},{0,2}};
					bondList[2]=new int [][]{{0,1}};
					bondList[4]=new int [][]{{2,3}};
					bondList[1]=new int [][]{{1,2}};
				}	else if (diagramIndex == 4){
					refBondList[0] = new int [][]{{0,1}};
					refBondList[1] = new int [][]{{1,2},{2,3}};
					bondList[5]=new int [][]{{3,0},{1,3},{0,2}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2},{2,3}};
				}	else if (diagramIndex == 6){
					refBondList[1] = new int [][]{{0,1}};
					refBondList[0] = new int [][]{{1,2},{2,3}};
					bondList[5]=new int [][]{{3,0},{1,3},{0,2}};
					bondList[2]=new int [][]{{0,1}};
					bondList[1]=new int [][]{{1,2},{2,3}};
				}
		}	else if (numDiagram == 28) {
			bondList[5]=new int [][]{{0,1},{1,2},{1,3},{3,0},{2,0},{2,3}};
			
		}	else if (numDiagram == 29) {
			bondList[5]=new int [][]{{1,2},{1,3},{3,0},{2,0},{2,3}};
			bondList[1]=new int [][]{{0,1}};
			
		}	else if (numDiagram == 30){
				if (diagramIndex == 1){
					bondList[5]=new int [][]{{2,3},{3,1},{2,0},{0,3}};
					bondList[1]=new int [][]{{0,1},{1,2}};
				}	else if (diagramIndex == 2){
					bondList[5]=new int [][]{{2,3},{3,1},{2,0},{0,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[4]=new int [][]{{1,2}};
				}
		}	else if (numDiagram == 34){
				if (diagramIndex == 1){
					refBondList[0] = new int [][]{{0,1},{1,2},{2,3}};
					bondList[5]=new int [][]{{3,0},{0,2},{1,3}};
					bondList[1]=new int [][]{{0,1},{1,2},{2,3}};
				}	else if (diagramIndex == 2){
					refBondList[0] = new int [][]{{0,1}};
					refBondList[1] = new int [][]{{2,3}};
					refBondList[2] = new int [][]{{1,2}};
					bondList[5]=new int [][]{{3,0},{0,2},{1,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{2,3}};
					bondList[4]=new int [][]{{1,2}};
				}
		}	else if (numDiagram == 36){
			bondList[5]=new int [][]{{1,2},{3,0},{0,2},{1,3}};
			bondList[1]=new int [][]{{0,1},{2,3}};
		}
			else {
			throw new RuntimeException("This is strange");
		}	
        System.out.println("B4HS: "+HSB[4]);
		ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody>3, eRef, true);
        if (numDiagram == 12||numDiagram ==19||numDiagram ==34){
			ClusterBonds refBonds = new ClusterBonds(4,refBondList, false);
			refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fCARef4mer, fACRef4mer, fBCRef4mer});
        }
		refCluster.setTemperature(temperature);
		clusters = (ClusterBonds[]) Arrays.addObject(clusters, new ClusterBonds(nBody, bondList, false));
		targetCluster = new ClusterSumPolarizableWertheimProduct4Pt(clusters, new double[]{1}, new MayerFunction[]{fR, FCA, FAC, FCB, FBC, eR, e});
		if (flip) {
			System.out.println("Flipping is applied");
			((ClusterSumPolarizableWertheimProduct4Pt) targetCluster).setDeltaCut(deltaCut);
			((ClusterSumPolarizableWertheimProduct4Pt) targetCluster).setCaching(false);
			targetCluster = new ClusterCoupledFlipped(targetCluster, space);
		}
		targetCluster.setTemperature(temperature);

		SpeciesWater4P species = new SpeciesWater4P(space);
		species.setConformation(new ConformationWaterGCPM(space));
		final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, temperature, refCluster, targetCluster);
		ConfigurationClusterWertheimGCPM4Pt configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(), pCA);
		if (numDiagram == 28) configuration.initializeCoordinatesER(sim.box[1]);
		if (numDiagram == 29) configuration.initializeCoordinates(sim.box[1]);
		if ((numDiagram == 16 || numDiagram == 30)) {
			if (diagramIndex == 1) {
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(), pCA, pCA);
			} else if (diagramIndex == 4) {
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(), pCA, pAC);
			} else if (diagramIndex == 3) {
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(), pAC, pCA);
			} else if (diagramIndex == 2){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pCA,pBC);
			}	
			configuration.initializeCoordinates2(sim.box[1]);
			}
		if (numDiagram == 12 || numDiagram == 19 ||numDiagram == 34){
			if (diagramIndex == 1){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pCA,pCA, pCA);
			}	else if (diagramIndex == 4){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pCA,pAC, pAC);
			}	else if (diagramIndex == 2){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pCA,pBC, pAC);
			}	else if (diagramIndex == 3){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pAC,pCA, pBC);
			}	else if (diagramIndex == 5){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pCA,pAC, pCA);
			}	else if (diagramIndex == 6){
				configuration = new ConfigurationClusterWertheimGCPM4Pt(space, sim.getRandom(),pAC,pCA, pCA);
			}
			configuration.initializeCoordinates3(sim.box[0]);	
			configuration.initializeCoordinates3(sim.box[1]);
			}
			
		sim.setAccumulatorBlockSize((int)numSteps*10);
		sim.integratorOS.setNumSubSteps(1000);	

		if (false) {
			sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
			sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
			((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
			((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
			((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
			((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);

			simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
			simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
			simGraphic.makeAndDisplayFrame();
			((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box[1]).canvas).setBackgroundColor(Color.WHITE);
			sim.integratorOS.setNumSubSteps(100000);
			sim.setAccumulatorBlockSize(1000);

			// if running interactively, set filename to null so that it doens't read
			// (or write) to a refpref file
			sim.initRefPref(null, 10, false);
			sim.equilibrate(null, 20);
			sim.getController().addActivity(new ActivityIntegrate2(sim.integratorOS));
			if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
				throw new RuntimeException("Oops");
			}

			return;
		}
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+numDiagram+"_"+diagramIndex+"_"+temperature+"_"+flip : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        //sim.setRefPref(50.0);
        sim.initRefPref(refFileName, numSteps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, numSteps/20);  
        
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
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
                System.out.println("abs average: "+ratioAndError[0]*HSB[nBody]+", error: "+ratioAndError[1]*Math.abs(HSB[nBody]));
			}
		};
		IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
		progressReportListener.setInterval((int) (numSteps / 10));
		sim.integratorOS.getEventManager().addListener(progressReportListener);

		sim.integratorOS.getMoveManager().setEquilibrating(false);
		sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), numSteps);

		System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
		System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

		sim.printResults(HSB[nBody]);
	}


	public static class WertheimAssociatingFluidParam extends ParameterBase {
		public double temperature = 630;//reduced temperature
		public double sigmaHSRef = 5.0;
		public long numSteps = 10000;
		public double associationEnergy = 2500.0;
		public int numDiagram = 16;
		public int diagramIndex = 3;
		public boolean flip = false;
		public boolean bondingAngleRestriction = true;
	}

}
