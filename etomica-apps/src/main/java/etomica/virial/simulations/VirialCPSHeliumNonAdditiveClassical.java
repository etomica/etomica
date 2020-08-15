/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.potential.P3CPSNonAdditiveHeSimplified;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSumNonAdditiveTrimerEnergy;
import etomica.virial.MayerFunctionSphericalThreeBody;
import etomica.virial.MayerFunctionThreeBody;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/* 
 * Adapted by Kate from VirialGCPM
 * 
 * Computes only the nonadditive component of either the third, fourth, or fifth virial coefficient for the
 * ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009). 
 * 
 * 
 */


public class VirialCPSHeliumNonAdditiveClassical {
	
	protected static Vector r0, r1, r2;
	protected static Vector r01Vec, r12Vec, r02Vec;
	
	public VirialCPSHeliumNonAdditiveClassical(Space space) {
        r0 = space.makeVector();
        r1 = space.makeVector();
        r2 = space.makeVector();
        r01Vec = space.makeVector();
        r12Vec = space.makeVector();
        r02Vec = space.makeVector();
	}

    public static void main(String[] args) {
    	
    	Space space = Space3D.getInstance();


        VirialParam params = new VirialParam();
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long blocks; int nullRegionMethod; int stepsPerBlock; long blocksEq; boolean adjustStepFreq; boolean simplifiedP3NonAdd; boolean total;
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperatureK = params.temperature;
            blocks = params.blocks;
            stepsPerBlock = params.stepsPerBlock;
            blocksEq = params.blocksEq;
            adjustStepFreq = params.adjustStepFreq;
            sigmaHSRef = params.sigmaHSRef;
            nullRegionMethod = params.nullRegionMethod;
            simplifiedP3NonAdd = params.simplifiedP3NonAdd;
            total = params.total;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 10) {
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
            simplifiedP3NonAdd = Boolean.parseBoolean(args[8]);
            total = Boolean.parseBoolean(args[9]);
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

        System.out.println();
        System.out.println("Reference: B"+nPoints+"HS = "+HSB[nPoints]+", sigmaHSRef = "+sigmaHSRef);
        if (total) {
        	System.out.println("Target: Helium-4 B"+nPoints+" at T="+temperatureK+ " K");
        } else {
        	System.out.println("Target: Helium-4 B"+nPoints+"NonAdd at T="+temperatureK+ " K");
        }
        System.out.println("null region method on nonadditive potential = "+nullRegionMethod);
        if (simplifiedP3NonAdd) {
        	System.out.println("Simplified nonadditive potential used");
        }
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        
        
        P2HePCKLJS p2 = new P2HePCKLJS(space);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2);
        
        
        final P3CPSNonAdditiveHe p3NonAdd = new P3CPSNonAdditiveHe(space);
        final P3CPSNonAdditiveHeSimplified p3NonAddSimplified = new P3CPSNonAdditiveHeSimplified(space);
        p3NonAdd.setNullRegionMethod(nullRegionMethod);
        p3NonAddSimplified.setNullRegionMethod(nullRegionMethod);
        
        MayerFunctionThreeBody f3Target = new MayerFunctionSphericalThreeBody(p3NonAdd);
        MayerFunctionThreeBody f3TargetSimplified = new MayerFunctionSphericalThreeBody(p3NonAddSimplified);
        
        VirialDiagrams diagrams = new VirialDiagrams(nPoints, true, false);
        diagrams.setDoReeHoover(true);
        diagrams.setDoShortcut(true);
        
    	
    	ClusterSumNonAdditiveTrimerEnergy targetCluster;
    	if (simplifiedP3NonAdd) {
    		targetCluster = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, p3NonAddSimplified, nPoints>3, false);
    	} else {
    		targetCluster = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, p3NonAdd, nPoints>3, false);
    	}
    	targetCluster.setNo72B2B3NonAdd(false);
    	targetCluster.setTemperature(temperature);
    	
    	ClusterAbstract targetClusterT;
    	if (simplifiedP3NonAdd) {
    		targetClusterT = diagrams.makeVirialCluster(fTarget, f3TargetSimplified, true);
    	} else {
    		targetClusterT = diagrams.makeVirialCluster(fTarget, f3Target, true);
    	}
    	targetClusterT.setTemperature(temperature);
    	
    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        refCluster.setTemperature(temperature);

        
        SpeciesSpheresMono species = new SpeciesSpheresMono(space, new ElementSimple("A"));
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, temperature, refCluster,total ? targetClusterT : targetCluster, false);
                
        
        sim.integratorOS.setAdjustStepFraction(adjustStepFreq);
        System.out.println("adjustStepFreq = " + adjustStepFreq);
        
        System.out.println(blocks*stepsPerBlock+" steps ("+blocks+" blocks of "+stepsPerBlock+" steps)");
        int numSubSteps = 1000; //steps per "overlap" block
        sim.integratorOS.setNumSubSteps(numSubSteps);
        System.out.println(numSubSteps+" steps per overlap-sampling block");
        long steps = stepsPerBlock*blocks/numSubSteps;
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
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(stepsPerBlock);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0, false);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
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
        
        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), steps);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 50.0;   // Kelvin
        public long blocks = 1000;  //NOT overlap blocks
        public int stepsPerBlock = 1000;
        public long blocksEq=1000; //NOT overlap steps
        public boolean adjustStepFreq = false;
        public double sigmaHSRef = 3;
        private int nullRegionMethod = 2; // What we have been using so far.
        public boolean writeRefPref;
        public boolean simplifiedP3NonAdd = false;
        public boolean total = true; //compute total virial coefficients
    }
}
