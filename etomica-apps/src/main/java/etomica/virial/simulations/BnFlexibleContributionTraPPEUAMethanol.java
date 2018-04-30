/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.box.Box;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.traPPE.MethanolPotentialHelper;
import etomica.models.traPPE.SpeciesMethanol;
import etomica.potential.P3BondAngle;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Computes (part of) the flexible contribution to the virial coefficient
 * for the TraPPE-UA model of methanol (http://www.chem.umn.edu/groups/siepmann/trappe/intro.php)
 *   
 * K.R.S. Shaul 2010
 *  
 */

public class BnFlexibleContributionTraPPEUAMethanol {

    // to control whether or not graphics are used:
    protected static boolean graphics = false;

    public static void main(String[] args) {

        VirialParam params = new VirialParam();

        if (args.length == 4 ) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	params.numMolecules = Integer.parseInt(args[0]);
        	params.temperature = Integer.parseInt(args[1]);
        	params.numSteps = Integer.parseInt(args[2]);
        	params.sigmaHSRef = Double.parseDouble(args[3]);

        } else if (args.length != 0){
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }

        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;

        // number of overlap sampling steps
        // for each overlap sampling step, the simulation boxes are allotted
        // 1000 attempts for MC moves, total
        long steps = params.numSteps;

        double sigmaHSRef = params.sigmaHSRef;

        if (numMolecules == 3 ) {
        	System.out.println("Overlap-sampling calculation for T1 (see Caracciolo et al. 2006)");
        } else if (numMolecules == 4 ) {
        	System.out.println("Overlap-sampling calculation for B4 + 1/8*I4 - 5/2*I2*T1 (see Caracciolo et al. 2006)");
        }

        System.out.println("");
        System.out.print("Target system: TraPPE-UA model of methanol");

        System.out.println("");
        System.out.println("");

        //****************************************************************************
        // Hard-sphere reference value
        // ****************************************************************************

        double B2HS = Standard.B2HS(sigmaHSRef);
        double B3HS = Standard.B3HS(sigmaHSRef);

        double ref;

        if (numMolecules == 3) {

	        ref = 4*B2HS*B2HS;

            System.out.println("Reference diagram: 4*B2*B2 for hard spheres = " + ref+ " Angstroms^6");

        } else if (numMolecules == 4) {

            ref = 16*B2HS*B2HS*B2HS - 9*B2HS*B3HS;

            System.out.println("Reference diagram: 16*B2HS*B2HS*B2HS - 9*B2HS*B3HS = " + ref + " Angstroms^9");

        } else {
        	throw new IllegalArgumentException("Cannot yet compute diagrams for that order of virial coefficient");
        }

        System.out.println("HS diameter: "+sigmaHSRef + " Angstroms");

        //****************************************************************************
        // Directives for overlap sampling
        //****************************************************************************

        System.out.println("");
        System.out.println("Temperature: "+temperature+" K");
        System.out.println("");

        temperature = Kelvin.UNIT.toSim(temperature);

        Space space = Space3D.getInstance();

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        // U_a_b is a pairwise potential (2 molecules, a and b, are involved).
        // The directives for calculation of U_a_b are provided later.
        PotentialGroup U_a_b = new PotentialGroup(2);

        MayerGeneral fTarget = new MayerGeneral(U_a_b);


        //****************************************************************************
        // Define target and reference clusters
        //****************************************************************************

        ClusterAbstract targetCluster;
        ClusterSum refCluster;

        if (numMolecules == 3) {

            //*********************************************
        	// T1 = t1a - t1b (Same T1 as Caracciolo et al. 2006)
        	//*********************************************

            ClusterBonds t1a   = new ClusterBonds(4, new int[][][]{{{0,1},{0,2}}});
	        ClusterBonds t1b   = new ClusterBonds(4, new int[][][]{{{0,1},{3,2}}});

            //Swapping zeroth and third molecules (to improve sampling efficiency)
	        ClusterBonds t1a03 = new ClusterBonds(4, new int[][][]{{{3,1},{3,2}}});
	        ClusterBonds t1b03 = new ClusterBonds(4, new int[][][]{{{3,1},{0,2}}});

            targetCluster = new ClusterSum(new ClusterBonds[]{t1a, t1b, t1a03, t1b03}, new double[]{0.5,-0.5,0.5, -0.5}, new MayerFunction[]{fTarget});

            refCluster = new ClusterSum(new ClusterBonds[]{t1a}, new double[]{1}, new MayerFunction[]{fRef});

        } else if (numMolecules == 4) {


            //*********************************************
        	// S2 = s2a - s2b (Similar to T2 of Caracciolo et al. 2006)
        	//*********************************************

            ClusterBonds s2a     = new ClusterBonds(5, new int[][][]{{{0,1},{0,2},{0,3}}});
        	ClusterBonds s2a04   = new ClusterBonds(5, new int[][][]{{{4,1},{4,2},{4,3}}});

            ClusterBonds s2b     = new ClusterBonds(5, new int[][][]{{{0,1},{0,2},{4,3}}});
            ClusterBonds s2b04   = new ClusterBonds(5, new int[][][]{{{4,1},{4,2},{0,3}}});

            ClusterBonds s2b13   = new ClusterBonds(5, new int[][][]{{{0,3},{0,2},{4,1}}});
            ClusterBonds s2b0413 = new ClusterBonds(5, new int[][][]{{{4,3},{4,2},{0,1}}});

            ClusterBonds s2b23   = new ClusterBonds(5, new int[][][]{{{0,1},{0,3},{4,2}}});
            ClusterBonds s2b0423 = new ClusterBonds(5, new int[][][]{{{4,1},{4,3},{0,2}}});

            //*********************************************
        	// S3 = s3a - s3b (Similar to T3 of Caracciolo et al. 2006)
        	//*********************************************

            ClusterBonds s3a   = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{0,3}}});
            ClusterBonds s3a04 = new ClusterBonds(5, new int[][][]{{{4,1},{1,2},{4,3}}});

            ClusterBonds s3b   = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{4,3}}});
            ClusterBonds s3b04 = new ClusterBonds(5, new int[][][]{{{4,1},{1,2},{0,3}}});

            //*********************************************
        	// T4 = t4a - t4b (Same T4 as Caracciolo et al. 2006)
        	//*********************************************

            ClusterBonds t4a   = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{0,2},{0,3}}});
        	ClusterBonds t4a04 = new ClusterBonds(5, new int[][][]{{{4,1},{1,2},{4,2},{4,3}}});

            ClusterBonds t4b   = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{0,2},{4,3}}});
            ClusterBonds t4b04 = new ClusterBonds(5, new int[][][]{{{4,1},{1,2},{4,2},{0,3}}});

            ClusterBonds[] diagrams = new ClusterBonds[]{s2a,s2a04,s2b,s2b04,s2b13,s2b0413,s2b23,s2b0423,s3a,s3a04,s3b,s3b04,t4a,t4a04,t4b,t4b04};

            double w = 0.75;
            double[] weights = new double[]{w/3.0,w/3.0,-w/9.0,-w/9.0,-w/9.0,-w/9.0,-w/9.0,-w/9.0, w,w,-w,-w, w,w,-w,-w};

            targetCluster = new ClusterSum(diagrams, weights, new MayerFunction[]{fTarget});

            refCluster = new ClusterSum(new ClusterBonds[]{s2a,s3a,t4a}, new double[]{0.5,1.5,1.5}, new MayerFunction[]{fRef});

        } else {
        	throw new IllegalArgumentException("Cannot yet compute diagrams for that order of virial coefficient");
        }


        System.out.println("ClusterCoupledFlipped NOT employed.");
        System.out.println("Permutations exchanging the root molecule (0) and its doppelganger (4) are employed.");
        System.out.println("For S2b, permutations exchanging 3 with 1 and 2 are employed.");
        System.out.println("");

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println(steps*1000+" total attempted MC moves, ("+steps+" blocks of 1000)");


        //PotentialMaster potentialMaster = new PotentialMaster(space);


        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2 (space,new SpeciesMethanol(space),
                temperature,refCluster,targetCluster, true); //use first constructor; no need for intramolecular movement MC trial
    	//sim.setRandom(new RandomNumberGenerator(8));
    	SpeciesMethanol species = (SpeciesMethanol)sim.getSpecies(0);
    	MethanolPotentialHelper.initPotential(space, species, U_a_b);
    	//potentialMaster.addPotential(U_a_b, new ISpecies[] {species,species} );


        // INTRAmolecular harmonic bending potential
        double thetaEq = 108.5*Math.PI/180;
        double kTheta = Kelvin.UNIT.toSim(55400); // force constant [=] K;
        PotentialGroup U_bend = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);

        P3BondAngle uBending = new P3BondAngle(space);

        uBending.setAngle(thetaEq);
        uBending.setEpsilon(kTheta);

        U_bend.addPotential(uBending, new Atomset3IteratorIndexList(new int[][] {{0,1,2}}));
        // integrators share a common potentialMaster.  so just add to one
        sim.integrators[1].getPotentialMaster().addPotential(U_bend,new ISpecies[]{species});


        //***********************************************************************************
        // Specify the doppelgaenger molecule(s) (4 copies 0 w/r to position and orientation)
        //***********************************************************************************

        if (numMolecules == 3) {

            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(new int[]{0,1,2,0});

            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(new int[]{0,1,2,0});

            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(new
            int[]{0,1,2,0});

            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(new
            int[]{0,1,2,0});

        } else {

            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(new int[]{0,1,2,3,0});

            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(new int[]{0,1,2,3,0});

            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(new int[]{0,1,2,3,0});

            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(new int[]{0,1,2,3,0});

        }


//         sim.integratorOS.setAdjustStepFreq(false);
//         sim.integratorOS.setStepFreq0(1);

        Box referenceBox = sim.box[0];
        Box targetBox = sim.box[1];

        sim.integratorOS.setNumSubSteps(1000); // Is this necessary?


        //****************************************************************************
	    // Make sure we find a good starting configuration
	    //****************************************************************************


        for (int j=0; j<1000 &&
	    	sim.box[1].getSampleCluster().value(sim.box[1]) == 0; j++) {
	    	sim.integrators[1].doStep();
	    }
	    if (sim.box[1].getSampleCluster().value(sim.box[1]) == 0) {
	    	throw new RuntimeException("could not find a configuration for target system");
	    }
	    sim.accumulators[1].reset();


        //****************************************************************************
        // Graphics
        //   true to run graphics (and not collect data)
        //   false to not run graphics (and collect data)
        //****************************************************************************


        if (graphics) {

            referenceBox.getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            targetBox.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(referenceBox).setShowBoundary(false);
            simGraphic.getDisplayBox(targetBox).setShowBoundary(false);

            simGraphic.getDisplayBox(referenceBox).setLabel("Reference-System Sampling");
            simGraphic.getDisplayBox(targetBox).setLabel("Target-System Sampling");

            // Create instances of ColorSchemeByType for reference and target simulations
            ColorSchemeByType colorScheme0 = (ColorSchemeByType) simGraphic.getDisplayBox(referenceBox).getColorScheme();
            ColorSchemeByType colorScheme1 = (ColorSchemeByType) simGraphic.getDisplayBox(targetBox).getColorScheme();


            // Create instances of the types of molecular sites
            AtomType typeCH3 = species.getCH3Type();
            AtomType typeO = species.getOType();
            AtomType typeH = species.getHType();

            // Set color of each site type for each simulation

            colorScheme0.setColor(typeCH3, Color.GRAY);
            colorScheme0.setColor(typeO, Color.RED);
            colorScheme0.setColor(typeH, Color.WHITE);

            colorScheme1.setColor(typeCH3, Color.GRAY);
            colorScheme1.setColor(typeO, Color.RED);
            colorScheme1.setColor(typeH, Color.WHITE);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 100);
                    sim.equilibrate(null, 200);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }


        //****************************************************************************
        // Let's get things started...
        //****************************************************************************

        // If running interactively, don't use a refPref file.
        String refFileName = args.length > 0 ? "refpref"+numMolecules+"_"+temperature : null;
        // This will either read the refPref in from a file or run a short simulation to find it.
        sim.initRefPref(refFileName, steps/40);
        // Run another short simulation to find MC move step sizes and maybe narrow in more on the best "refPref."
        // If it does continue looking for a refPref, it will write the value to the file.
        sim.equilibrate(refFileName, steps/20);

        sim.setAccumulatorBlockSize(steps);

        System.out.println();
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        System.out.println();

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());


        sim.printResults(ref);

    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

        // number of molecules in simulation (e.g., 2 for B2 calculation)
    	public int numMolecules = 5;

        public double temperature = 300.0;   // Kelvin

        // number of overlap sampling steps
        // for each overlap sampling step, the simulation boxes are allotted
        // 1000 attempts for MC moves, total
        public long numSteps = 1000;

        public double sigmaHSRef = 8.0;

    }
    
}


