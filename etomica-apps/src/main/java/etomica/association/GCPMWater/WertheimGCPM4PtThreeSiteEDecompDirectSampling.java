/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association.GCPMWater;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.P2HardAssociationGCPMReference;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.models.water.SpeciesWater4P;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.simulations.SimulationVirial;

import java.awt.*;

/**
 * repulsive potential: energy of pair is greater than -association E(cal/mol)
 * attractive potential: energy of pair is less than -association E(cal/mol)
 * Wertheim's three attraction-site model of GCPM
 * Direct Sampling
 * Non-additive 4th order diagram
 *
 * @author Hye Min Kim
 */

public class WertheimGCPM4PtThreeSiteEDecompDirectSampling {

	public static void main(String[] args) {
		VirialAssociatingFluidParam params = new VirialAssociatingFluidParam();
		Space space = Space3D.getInstance();
		ParseArgs.doParseArgs(params, args);
		
		final int nBody = 4;//W4
		int numDiagram = params.numDiagram;
		int diagramIndex = params.diagramIndex;
		double temperature = params.temperature;
		double associationEnergy = params.associationEnergy;
		long numSteps = params.numSteps;
		boolean bondingAngleRestriction = params.bondingAngleRestriction;

		final double[] HSB = new double[9];
        HSB[2] = -35.238;//hard-chain reference
        HSB[3] = -35.238*-35.238;//hard-chain reference
        HSB[4] = -35.238*-35.238*-35.238;//hard-chain reference
        System.out.println("Direct sampling");
        System.out.println("Wertheim's non-additive Wertheim diagram, final version, fij in delta B4 is NOT decomposed into fR and fAssociation");
        System.out.println("Wertheim's three association site model (2 Hydrogens and 1 Oxygen)");
        System.out.println("repulsive potential: greater than -" +associationEnergy+ " cal/mol, attractive potential: less than -"+associationEnergy+" cal/mol");
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("Diagram "+numDiagram+" "+diagramIndex+"E "+numSteps+" numSteps");
        System.out.println("temperature: "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        System.out.println("association Energy: "+associationEnergy+"cal/mole");
        if (bondingAngleRestriction){
        	System.out.println("This calculation includes bonding angle criteria(cosTheta1<-0.5,cosTheta3<0)");
        }
    	Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	associationEnergy = calPerMole.toSim(associationEnergy);
    	
        IPotentialMolecular pTarget = new PNWaterGCPM(space);
        
        MayerEGeneral e = new MayerEGeneral(pTarget);

		PNWaterGCPMThreeSite pR = new PNWaterGCPMThreeSite(space,associationEnergy,false,bondingAngleRestriction); //GCPM potential
		PNWaterGCPMThreeSite pAC = new PNWaterGCPMThreeSite(space,associationEnergy,true,bondingAngleRestriction); //assciation between site A on mol1 and C on mol2
		pAC.setBondType(1);
		PNWaterGCPMThreeSite pBC = new PNWaterGCPMThreeSite(space,associationEnergy,true,bondingAngleRestriction); //assciation between site B on mol1 and C on mol2
		pBC.setBondType(2);
		PNWaterGCPMThreeSite pCA = new PNWaterGCPMThreeSite(space,associationEnergy,true,bondingAngleRestriction); //assciation between site C on mol1 and A on mol2
		pCA.setBondType(3);
		PNWaterGCPMThreeSite pCB = new PNWaterGCPMThreeSite(space,associationEnergy,true,bondingAngleRestriction); //assciation between site C on mol1 and A on mol2
		pCB.setBondType(4);
		MayerEGeneral eR = new MayerEGeneral(pR);//repulsion eR function
		MayerGeneral fR = new MayerGeneral(pR);//Mayer f function of reference part
		MayerGeneral fAC = new MayerGeneral(pAC);//Mayer f function of association part
		MayerGeneral fBC = new MayerGeneral(pBC);//Mayer f function of association part
		MayerGeneral fCA = new MayerGeneral(pCA);//Mayer f function of association part
		MayerGeneral fCB = new MayerGeneral(pCB);//Mayer f function of association part
		MayerFunctionProductGeneral FAC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fAC}, new double[]{1,1});//association bond between site A on mol1 and C on mol2
		MayerFunctionProductGeneral FBC = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fBC}, new double[]{1,1});//association bond between site B on mol1 and C on mol2
		MayerFunctionProductGeneral FCA = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fCA}, new double[]{1,1});//association bond between site C on mol1 and A on mol2
		MayerFunctionProductGeneral FCB = new MayerFunctionProductGeneral(space, new MayerFunction[]{eR,fCB}, new double[]{1,1});//association bond between site C on mol1 and A on mol2
        P2HardAssociationGCPMReference pCARef = new P2HardAssociationGCPMReference(space, true);
        pCARef.setBondType(3);
		MayerGeneral fCARef = new MayerGeneral(pCARef);
        MayerEGeneral eCARef = new MayerEGeneral(pCARef);
        P2HardAssociationGCPMReference pACRef = new P2HardAssociationGCPMReference(space, true);
        pACRef.setBondType(1);
		MayerGeneral fACRef = new MayerGeneral(pACRef);
        MayerEGeneral eACRef = new MayerEGeneral(pACRef);
        P2HardAssociationGCPMReference pBCRef = new P2HardAssociationGCPMReference(space, true);
        pBCRef.setBondType(2);
		MayerGeneral fBCRef = new MayerGeneral(pBCRef);
        MayerEGeneral eBCRef = new MayerEGeneral(pBCRef);
		
		int nBondTypes = 6;//fR,FCA,FAC,FCB,,FBC,eR
		ClusterBonds[] clusters = new ClusterBonds[0];
		int[][][] bondList = new int[nBondTypes][][];
		int [][][]refBondList = new int[nBondTypes][][];
		ClusterAbstract targetCluster = null;
		if (numDiagram == 12){			
				if (diagramIndex == 5){
					refBondList[0] = new int [][]{{0,1},{2,3}};
					refBondList[1] = new int [][]{{1,2}};
					bondList[0]=new int [][]{{3,0}};
					bondList[1]=new int [][]{{0,1},{2,3}};
					bondList[2]=new int [][]{{1,2}};
					bondList[5]=new int [][]{{0,2},{1,3}};
				}
		}	else if (numDiagram == 19){
				if (diagramIndex == 3){
					refBondList[1] = new int [][]{{0,1},{2,3}};
					refBondList[0] = new int [][]{{1,2}};		
					bondList[0]=new int [][]{{3,0},{1,3}};
					bondList[2]=new int [][]{{0,1},{2,3}};
					bondList[1]=new int [][]{{1,2}};
					bondList[3]=new int [][]{{0,2}};
				}	else if (diagramIndex == 4){
					refBondList[0] = new int [][]{{0,1}};
					refBondList[1] = new int [][]{{1,2},{2,3}};		
					bondList[0]=new int [][]{{3,0},{1,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2},{2,3}};
					bondList[3]=new int [][]{{0,2}};
				}	else if (diagramIndex == 6){
					refBondList[1] = new int [][]{{0,1}};
					refBondList[0] = new int [][]{{1,2},{2,3}};
					bondList[0]=new int [][]{{3,0},{1,3}};
					bondList[2]=new int [][]{{0,1}};
					bondList[1]=new int [][]{{1,2},{2,3}};
					bondList[3]=new int [][]{{0,2}};
				}
		}	else if (numDiagram == 23){
				if (diagramIndex == 2){
					refBondList[0] = new int [][]{{1,2},{1,3}};
					refBondList[1] = new int [][]{{0,1}};
					bondList[5]=new int [][]{{2,3},{0,2},{0,3}};
					bondList[1]=new int [][]{{1,2},{1,3}};
					bondList[2]=new int [][]{{0,1}};
				}	else if (diagramIndex == 4){
					refBondList[0] = new int [][]{{0,1},{1,3}};
					refBondList[1] = new int [][]{{1,2}};
					bondList[5]=new int [][]{{2,3},{0,3},{0,2}};
					bondList[1]=new int [][]{{0,1},{1,3}};
					bondList[2]=new int [][]{{1,2}};
				}	else if (diagramIndex == 5){
					refBondList[0] = new int [][]{{0,1}};
					refBondList[1] = new int [][]{{1,2},{1,3}};
					bondList[5]=new int [][]{{2,3},{0,3},{0,2}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2},{1,3}};
				}			
		}	else if (numDiagram == 30){
				if (diagramIndex == 1){
					bondList[0]=new int [][]{{2,3},{3,1},{2,0},{0,3}};
					bondList[1]=new int [][]{{0,1},{1,2}};
				}	else if (diagramIndex == 2){
					bondList[0]=new int [][]{{2,3},{3,1},{2,0},{0,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2}};
				}	else if (diagramIndex == 3){
					bondList[0]=new int [][]{{2,3},{3,1},{0,3}};
					bondList[2]=new int [][]{{0,1}};
					bondList[1]=new int [][]{{1,2}};
					bondList[3]=new int [][]{{0,2}};
				}	else if (diagramIndex == 4){
					bondList[0]=new int [][]{{2,3},{3,1},{0,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{1,2}};
					bondList[3]=new int [][]{{0,2}};
				}
		}	else if (numDiagram == 32){
			if (diagramIndex == 1){
				refBondList[0] = new int [][]{{0,1},{1,3}};
				refBondList[2] = new int [][]{{1,2}};
				bondList[5]=new int [][]{{2,3},{3,0},{0,2}};
				bondList[1]=new int [][]{{0,1},{1,3}};
				bondList[4]=new int [][]{{1,2}};
			}	else if (diagramIndex == 3){
				refBondList[0] = new int [][]{{1,3}};
				refBondList[1] = new int [][]{{0,1},{1,2}};
				bondList[5]=new int [][]{{2,3},{0,2},{0,3}};
				bondList[1]=new int [][]{{1,3}};
				bondList[2]=new int [][]{{0,1},{1,2}};
			}	else if (diagramIndex == 6){
				refBondList[0] = new int [][]{{0,1}};
				refBondList[1] = new int [][]{{1,3}};
				refBondList[2] = new int [][]{{1,2}};
				bondList[5]=new int [][]{{2,3},{0,2},{0,3}};
				bondList[1]=new int [][]{{0,1}};
				bondList[2]=new int [][]{{1,3}};
				bondList[4]=new int [][]{{1,2}};
			}
	}	else if (numDiagram == 34){
				if (diagramIndex == 1){
					refBondList[0] = new int [][]{{0,1},{1,2},{2,3}};	
					bondList[0]=new int [][]{{3,0},{0,2},{1,3}};
					bondList[1]=new int [][]{{0,1},{1,2},{2,3}};
				}	else if (diagramIndex == 2){
					refBondList[0] = new int [][]{{0,1}};
					refBondList[1] = new int [][]{{2,3}};
					refBondList[2] = new int [][]{{1,2}};
					bondList[0]=new int [][]{{3,0},{0,2},{1,3}};
					bondList[1]=new int [][]{{0,1}};
					bondList[2]=new int [][]{{2,3}};
					bondList[4]=new int [][]{{1,2}};
				}		
		}	else if (numDiagram == 36){
			bondList[0]=new int [][]{{1,2},{3,0},{0,2},{1,3}};
			bondList[1]=new int [][]{{0,1},{2,3}};
		}
			else {
			throw new RuntimeException("This is strange");
		}	
		
		ClusterBonds refBonds = new ClusterBonds(nBody,refBondList, false);
		ClusterAbstract refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double []{1}, new MayerFunction[]{fCARef,fACRef,fBCRef});
		refCluster.setTemperature(temperature);
		clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
		targetCluster = new ClusterSumPolarizableWertheimProduct4Pt(clusters,new double []{1}, new MayerFunction[]{fR,FCA,FAC,FCB,FBC,eR,e});
        targetCluster.setTemperature(temperature);
		final SimulationVirial sim = new SimulationVirial(space, new SpeciesFactoryWaterGCPM(), temperature,ClusterWeightAbs.makeWeightCluster(refCluster), refCluster, new ClusterAbstract[]{targetCluster});
		ConfigurationClusterWertheimGCPMDirectSampling	configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pCA, pCA);
		if ((numDiagram == 12 || numDiagram == 34) && diagramIndex == 2){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pBC, pAC);
		}
		if (numDiagram == 34 && diagramIndex == 3){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pAC, pCA, pBC);
		}
		if ((numDiagram == 13 && diagramIndex == 3)||(numDiagram == 20 && diagramIndex == 4)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pCA, pAC);
		}
		if ((numDiagram == 13 && diagramIndex == 4)||(numDiagram == 20 && diagramIndex == 6)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pBC, pCA);
		}
		if ((numDiagram == 13 && diagramIndex == 5)||(numDiagram == 20 && diagramIndex == 8)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pBC, pAC);
		}
		if ((numDiagram == 13 && diagramIndex == 7)||(numDiagram == 20 && diagramIndex == 11)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pCA, pBC);
		}
		if ((numDiagram == 13 && (diagramIndex == 8||diagramIndex == 9))||((numDiagram == 12||numDiagram == 26) && diagramIndex == 5)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pAC, pCA);
		}
		if ((numDiagram == 13 && diagramIndex == 10)||(numDiagram == 26 && diagramIndex == 4)){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pAC, pAC);
		}
		if (numDiagram == 26 && diagramIndex == 6){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pAC, pCA, pCA);
		}
		configuration.initializeCoordinatesTetramer(sim.box);
		if ((numDiagram == 23 ||numDiagram == 32) && diagramIndex == 1){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pBC, pCA);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		if (numDiagram == 23 && diagramIndex == 2){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pAC, pCA, pCA);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		if ((numDiagram == 23 ||numDiagram == 32) && diagramIndex == 3){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pAC, pAC, pCA);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		if (numDiagram == 23 && diagramIndex == 4){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pAC, pCA);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		if (numDiagram == 23 && diagramIndex == 5){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pAC, pAC);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		if ((numDiagram == 23 ||numDiagram == 32) && diagramIndex == 6){	
			configuration = new ConfigurationClusterWertheimGCPMDirectSampling(space, sim.getRandom(), pCA, pBC, pAC);
			configuration.initializeCoordinatesBranch(sim.box);
		}
		



		sim.accumulator.setBlockSize(numSteps/100);
	

		if (true) {
            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ISpecies species = sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(1), Color.RED);

            simGraphic.getDisplayBox(sim.box).setShowBoundary(false);
            simGraphic.makeAndDisplayFrame();
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            return;
        }
		sim.equilibrate(numSteps/40);
ActivityIntegrate ai = new ActivityIntegrate(sim.integrator, numSteps);
System.out.println("equilibration finished");
        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
sim.getController().runActivityBlocking(ai);
        
        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        double referenceAverage = ((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0];
        double referenceError = ((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0];
        double targetAverage = ((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1];
        double targetError = ((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1];
        double ratio = ((DataGroup)sim.accumulator.getData()).getData(sim.accumulator.RATIO.index).getValue(1);
        double error = ((DataGroup)sim.accumulator.getData()).getData(sim.accumulator.RATIO_ERROR.index).getValue(1);
        System.out.println("reference average: "+referenceAverage+", error: "+referenceError);
		System.out.println("target average: "+targetAverage+", error: "+targetError);
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("diagram value: "+ratio*HSB[nBody]+", error: "+Math.abs(error*HSB[nBody]));
    }
    	
	public static class VirialAssociatingFluidParam extends ParameterBase {
		public double temperature = 630;//Kelvin
		public long numSteps = 100000;
		public double associationEnergy = 2500.0;
		public int numDiagram = 23;
		public int diagramIndex = 2;
		public boolean bondingAngleRestriction = false;
		
	}

}
