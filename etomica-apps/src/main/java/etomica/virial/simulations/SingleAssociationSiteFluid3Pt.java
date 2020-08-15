/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2HardAssociationCone;
import etomica.potential.P2LennardJones;
import etomica.potential.P2MoleculeMonatomic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * repulsive potential: Lennard-Jones potential
 * attractive potential: short-handed square well potential 
 * Wertheim's single attraction-site model
 * 2points and 3 points diagram calculation
 *
 * @author Hye Min Kim
 */

public class SingleAssociationSiteFluid3Pt {

	public static void main(String[] args) {
		VirialAssociatingFluidParam params = new VirialAssociatingFluidParam();
		Space space = Space3D.getInstance();
		
		int rhopoint=params.rhopoint;
		int rho0point=params.rho0point;
		int diagramIndex = params.diagramIndex;
		final int nBody = rhopoint + rho0point;
		double temperature=params.temperature;
		double sigmaHSRef = params.sigmaHSRef;
		double wellConstant = params.wellConstant;
		long numSteps = params.numSteps;
		
		if (args.length == 7) {
        	temperature = Double.parseDouble(args[0]);
            sigmaHSRef = Double.parseDouble(args[1]);
            numSteps = Integer.parseInt(args[2]);
            wellConstant = Double.parseDouble(args[3]);
            rhopoint = Integer.parseInt(args[4]);
            rho0point = Integer.parseInt(args[5]);
            diagramIndex = Integer.parseInt(args[6]);
            
        } else if (args.length != 0){
        	throw new IllegalArgumentException("Wrong number of arguments");
        }
		double sigma=1;
		double epsilon=1;
		final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("Square-Well epsilon =" +wellConstant+ "*epsilon");
        System.out.println("repulsive potential : Lennard-Jones potential");
        System.out.println("attractive potential : short-handed square well potential");
		P2HardAssociationCone p = new P2HardAssociationCone(space, sigma, epsilon, Double.POSITIVE_INFINITY, wellConstant); //Lennard-Jones potential+square-well site-site attraction potential, 1.0=cutoffFactor
		P2MoleculeMonatomic pMolecule = new P2MoleculeMonatomic(p);
		P2LennardJones pR = new P2LennardJones(space, sigma, epsilon);//L-J potential
		MayerGeneralSpherical fR = new MayerGeneralSpherical(pR);//repulsion Mayer fR function
		MayerESpherical eR = new MayerESpherical(pR);//repulsion eR function
		MayerGeneral f = new MayerGeneral(pMolecule);//usual Mayer f function
		MayerFunctionSumGeneral F = new MayerFunctionSumGeneral(space, new MayerFunction[]{f,fR}, new double[]{1,-1});//F=f-fR
		
		MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
		
		int nBondTypes = 3;//fR,F,eR
		ClusterBonds[] clusters = new ClusterBonds[0];
		int[][][] bondList = new int[nBondTypes][][];	
        ClusterSum targetCluster = null;
            		
		  if (rhopoint == 3 && rho0point == 0) {
			bondList[0]=new int [][]{{0,1},{1,2},{2,0}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,F,eR});
			
		}  else if (rho0point == 3 && rhopoint == 0){
			if (diagramIndex == 0) {
			bondList[1] = new int [][]{{0,1},{1,2}};
			bondList[2] = new int [][]{{2,0}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,F,eR});
			}
			else if (diagramIndex == 1) {
				bondList[1] = new int [][]{{0,1},{1,2},{2,0}};
				clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
				targetCluster = new ClusterSum(clusters,new double []{1}, new MayerFunction[]{fR,F,eR});
			}
			
		}  else if (rho0point == 2 && rhopoint ==1){
			bondList[0] = new int [][]{{1,2},{2,0}};
			bondList[1] = new int [][]{{0,1}};
			clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
			  targetCluster = new ClusterSum(clusters, new double[]{1}, new MayerFunction[]{fR, F, eR});
		  } else {
			  throw new RuntimeException("This is strange");
		  }


		ClusterAbstract refCluster = Standard.virialCluster(nBody, fRef, nBody > 3, eRef, true);
		refCluster.setTemperature(temperature);
		targetCluster.setTemperature(temperature);
		final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresRotating(space, new ElementSimple("O")), temperature, refCluster, targetCluster);
		ConfigurationClusterMove configuration = new ConfigurationClusterMove(space, sim.getRandom());
		configuration.initializeCoordinates(sim.box[1]);
		sim.setAccumulatorBlockSize((int) numSteps * 10);
		sim.integratorOS.setNumSubSteps(1000);

		if (false) {
			sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
			sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
			simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
			simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
			Species species = (Species) sim.getSpecies(0);
			AtomType typeLJ = species.getAtomType(0);
			simGraphic.makeAndDisplayFrame();

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
        String refFileName = args.length > 0 ? "refpref"+rhopoint+"_"+rho0point+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, numSteps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, numSteps/40);
                
        System.out.println("equilibration finished");

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
		for (int i = 0; i < 2; i++) {
			System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
		}
sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), numSteps);

		System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());

		sim.printResults(HSB[nBody]);
	}
    
	
	public static class VirialAssociatingFluidParam extends ParameterBase {
		public double temperature = 1.5;//reduced temperature
		public double sigmaHSRef = 1.5;
		public long numSteps = 100000;
		public double wellConstant = 16.0;
		public int rhopoint = 3;
		public int rho0point = 0;
		public int diagramIndex = 0;
		
	}

}
