/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.potential.P3CPSNonAdditiveHeSimplified;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/*
 * Adapted by Kate from VirialCPSHeliumNonAdditiveClassical
 *
 * Computes only the nonadditive component of either the third, fourth, or fifth classical virial coefficient for the
 * ab initio non-additive trimer potential for He-4 developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009)
 * MINUS that for a simplified version of this potential.
 *
 *
 */


public class VirialCPSHeliumNonAdditiveClassical_Correction {
	
	protected static Vector r0, r1, r2;
	protected static Vector r01Vec, r12Vec, r02Vec;
	
	public VirialCPSHeliumNonAdditiveClassical_Correction(Space space) {
        r0 = space.makeVector();
        r1 = space.makeVector();
        r2 = space.makeVector();
        r01Vec = space.makeVector();
        r12Vec = space.makeVector();
        r02Vec = space.makeVector();
	}

    public static void main(String[] args) {
    	
    	Space space = Space3D.getInstance();

    	VirialCPSHeliumNonAdditiveClassical_Correction calc = new VirialCPSHeliumNonAdditiveClassical_Correction(space);
      

        VirialParam params = new VirialParam();
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long blocks; int nullRegionMethod; int stepsPerBlock; long blocksEq; boolean adjustStepFreq; 
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperatureK = params.temperature;
            blocks = params.blocks;
            stepsPerBlock = params.stepsPerBlock;
            blocksEq = params.blocksEq;
            adjustStepFreq = params.adjustStepFreq;
            sigmaHSRef = params.sigmaHSRef;
            nullRegionMethod = params.nullRegionMethod;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 8) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperatureK = Double.parseDouble(args[1]);
            blocks = Integer.parseInt(args[2]);
            stepsPerBlock = Integer.parseInt(args[3]);
            blocksEq = Integer.parseInt(args[4]);
            adjustStepFreq = Boolean.parseBoolean(args[5]);
            sigmaHSRef = Double.parseDouble(args[6]);
            nullRegionMethod = Integer.parseInt(args[7]);
            params.writeRefPref = true;
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }

        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Helium overlap sampling B"+nPoints+"NonAdd at T="+temperatureK+ " K");
        System.out.println("null region method = "+nullRegionMethod);
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        
        
        P2HePCKLJS p2 = new P2HePCKLJS(space);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2);
        
        final P3CPSNonAdditiveHe p3NonAdd = new P3CPSNonAdditiveHe(space);
        p3NonAdd.setNullRegionMethod(nullRegionMethod);
    	ClusterSumNonAdditiveTrimerEnergy targetCluster1 = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, p3NonAdd, nPoints>3, false);
    	targetCluster1.setNo72B2B3NonAdd(false);
    	targetCluster1.setTemperature(temperature);
    	
    	final P3CPSNonAdditiveHeSimplified p3NonAddSimplified = new P3CPSNonAdditiveHeSimplified(space);
        p3NonAddSimplified.setNullRegionMethod(nullRegionMethod);
        ClusterSumNonAdditiveTrimerEnergy targetCluster2 = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, p3NonAddSimplified, nPoints>3, false);
    	targetCluster2.setNo72B2B3NonAdd(false);
    	targetCluster2.setTemperature(temperature);
        
    	ClusterAbstract[] targetDiagrams2 = new ClusterAbstract[1];
    	targetDiagrams2[0] = targetCluster2;
    	ClusterDifference targetCluster = new ClusterDifference(targetCluster1, targetDiagrams2);
    	
    
    	
    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        refCluster.setTemperature(temperature);


        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")),
                temperature, refCluster, targetCluster, false);
        
        System.out.println(blocks*stepsPerBlock+" steps ("+blocks+" blocks of "+stepsPerBlock+" steps)");
        int numSubSteps = 1000; //steps per "overlap" block
        sim.integratorOS.setNumSubSteps(numSubSteps);
        System.out.println(numSubSteps+" steps per overlap-sampling block");
        long steps = stepsPerBlock*blocks/numSubSteps;

        sim.integratorOS.setAdjustStepFraction(adjustStepFreq);
        System.out.println("adjustStepFreq = " + adjustStepFreq);
        
        sim.setAccumulatorBlockSize(stepsPerBlock);
        
        ///////////////////////////////////////////////
        // Initialize non-overlapped configuration
        ///////////////////////////////////////////////
        
        IAtomList atoms = sim.box[1].getLeafList();
        if (nPoints == 3) {
	        for (int i = 1; i<atoms.size(); i++) {
	        	atoms.get(i).getPosition().setX(0, i*sigmaHSRef);
	        }
        } else if (nPoints == 4) {
	        
	        atoms.get(1).getPosition().setX(0, sigmaHSRef);
	        
	        atoms.get(2).getPosition().setX(0, sigmaHSRef);
	        atoms.get(2).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.get(3).getPosition().setX(1, sigmaHSRef);
	        
        } else if (nPoints == 5) {
        	
        	atoms.get(1).getPosition().setX(0, sigmaHSRef);
        	atoms.get(1).getPosition().setX(1, sigmaHSRef);
        	
	        atoms.get(2).getPosition().setX(0, sigmaHSRef);
	        atoms.get(2).getPosition().setX(1, -sigmaHSRef);
	        
	        atoms.get(3).getPosition().setX(0, -sigmaHSRef);
	        atoms.get(3).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.get(4).getPosition().setX(0, -sigmaHSRef);
	        atoms.get(4).getPosition().setX(1, -sigmaHSRef);
	        
        } else {
        	throw new RuntimeException("Wrong number of points");
        }
        
        /*
        IAtomList atoms0 = sim.box[0].getLeafList();
        for (int i=1;i<atoms0.getAtomCount();i++) {
        	atoms0.getAtom(i).getPosition().setX(0, i*sigmaHSRef*1.3);
        } */ 
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesSpheresMono species = (SpeciesSpheresMono)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0, false);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/blocks*blocksEq); // 5000 IntegratorOverlap steps = 5e6 steps
        System.out.println((stepsPerBlock*blocksEq) + " equilibration steps ("+blocksEq+" blocks of "+stepsPerBlock+" steps)"); 
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize(stepsPerBlock);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());
        
        IAction progressReport = new IAction() {
            public void actionPerformed() {
                //System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                IAtomList atoms = sim.box[1].getLeafList();
                
                r0.E( atoms.get(0).getPosition() );
                r1.E( atoms.get(1).getPosition() );
                r2.E( atoms.get(2).getPosition() );
                
                r01Vec.Ev1Mv2(r0,r1);
                r02Vec.Ev1Mv2(r0,r2);
                r12Vec.Ev1Mv2(r1,r2);
                
	    		double r01 = (Math.sqrt(r01Vec.squared()));
	    		double r02 = (Math.sqrt(r02Vec.squared()));
	    		double r12 = (Math.sqrt(r12Vec.squared()));
	    		 
	    		
	    		double U = Kelvin.UNIT.fromSim(p3NonAdd.energy(atoms));
	    		 System.out.println(r01+"  "+r02+"  " +r12+"  " +U);

            }
        };

        //sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(progressReport, 1 ));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
    }



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
    	public int nPoints = 3;
        public double temperature = 50.0;   // Kelvin
        public long blocks = 1000;  //NOT overlap blocks
        public int stepsPerBlock = 1000;
        public long blocksEq=1000; //NOT overlap steps
        public boolean adjustStepFreq = false;
        public double sigmaHSRef = 3;
        private int nullRegionMethod = 2; // What we have been using so far.
        public boolean writeRefPref;
        public boolean simplifiedP3NonAdd = true;
    }
}
