/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.listener.IntegratorListenerAction;
import etomica.models.traPPE.MethanolPotentialHelper;
import etomica.models.traPPE.SpeciesFactoryMethanol;
import etomica.models.traPPE.SpeciesMethanol;
import etomica.potential.P3BondAngle;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Mayer-sampling MC simulation for methanol or ethanol using the TraPPE model: http://www.chem.umn.edu/groups/siepmann/trappe/intro.php
 *   
 * Set species, model (with or without point charges), virial-coefficent type (B2, B3, etc.),
 * and temperature in VirialParam() method at bottom of file. 
 *   
 * K.R. Schadel 2008
 *  
 */



public class VirialTraPPEAlcohol {

    // ethanol = false: methanol
    // ethanol = true: ethanol
    protected static boolean ethanol = false;
    protected static boolean graphics = false;
    
    public static void main(String[] args) {

        VirialParam params = new VirialParam();

        /*System.out.println(""+ (Double.POSITIVE_INFINITY-Double.POSITIVE_INFINITY));
        System.out.println(""+ Math.sqrt(Double.POSITIVE_INFINITY));
        System.out.println("" + (Double.POSITIVE_INFINITY)*(Double.POSITIVE_INFINITY));
        System.exit(1);*/

        // enables one to overwrite parameters values in VirialParam() and use those provided in string instead
        if (args.length ==3 ) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	params.numMolecules = Integer.parseInt(args[0]);
        	params.temperature = Integer.parseInt(args[1]);
        	params.numSteps = Integer.parseInt(args[2]);
        } else if (args.length != 0){
        	throw new IllegalArgumentException("Incorrect number of arguments passed to VirialRowleyAlcohol.");
        }


        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;

        // number of overlap sampling steps
        // for each overlap sampling step, the simulation boxes are allotted
        // 1000 attempts for MC moves, total
        long steps = params.numSteps;


        // Diameter of hard spheres in reference system
        // Should be about the size of the molecules in the target system
        double sigmaHSRef;
        if (ethanol) {
        	sigmaHSRef = 5;
        }
        else {
            sigmaHSRef = 8;
        }

        final double[] HSB = new double[8];

        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);

        System.out.println();
        System.out.println("sigmaHSRef: "+sigmaHSRef + " Angstroms");
        System.out.println("B2HS: "+HSB[2] + " Angstroms^3");
        System.out.println("B3HS: "+HSB[3] + " Angstroms^6"); // + " = " +(HSB[3]/(HSB[2]*HSB[2]))+ " B2HS^2");
        System.out.println("B4HS: "+HSB[4] + " Angstroms^9"); // + " = " +(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5] + " Angstroms^12"); // What's with this?: + " = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6] + " Angstroms^15"); //  + " = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7] + " Angstroms^18"); // + " = 0.013046 B2HS^6");

        Space space = Space3D.getInstance();


        /*
        ****************************************************************************
        ****************************************************************************
        Directives for overlap sampling
        ****************************************************************************
        ****************************************************************************
        */

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);

        // U_a_b is a pairwise potential (2 molecules, a and b, are involved).
        // The directives for calculation of U_a_b are provided later.
        PotentialGroup U_a_b = new PotentialGroup(2);


        MayerGeneral fTarget = new MayerGeneral(U_a_b);
        MayerEGeneral eTarget = new MayerEGeneral(U_a_b);

        System.out.println();
        System.out.print("TraPPE model for");
        if (ethanol) {
        	System.out.print(" ethanol");
        }
        else {
        	System.out.print(" methanol");
        }

        System.out.println();
        System.out.println("B"+numMolecules+" at "+temperature+"K");

        temperature = Kelvin.UNIT.toSim(temperature); // What are the simulation units for T?

        ClusterAbstract targetCluster = Standard.virialCluster(numMolecules, fTarget, numMolecules>3, eTarget, true);

        // These models have point charges
        ((ClusterSum)targetCluster).setCaching(false);
        targetCluster = new ClusterCoupledFlipped(targetCluster, space);

        targetCluster.setTemperature(temperature);

        ClusterAbstract refCluster = Standard.virialCluster(numMolecules, fRef, numMolecules>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println(steps*1000+" total attempted MC moves ("+steps+" blocks of 1000)");


        //PotentialMaster potentialMaster = new PotentialMaster(space);

        /*if(ethanol) {
        	sim = new SimulationVirialOverlap (space,new SpeciesFactoryEthanol(),
        			temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesEthanol species = (SpeciesEthanol)sim.species;
        	EthanolPotentialHelper.initPotential(space, species, U_a_b, pointCharges);
        	potentialMaster.addPotential(U_a_b, new ISpecies[] {species,species} );
        }
        else {*/
        	final SimulationVirialOverlap sim;
        	sim = new SimulationVirialOverlap (space,new SpeciesFactoryMethanol(),
                    temperature,refCluster,targetCluster, true); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesMethanol species = (SpeciesMethanol)sim.getSpecies(0);
        	MethanolPotentialHelper.initPotential(space, species, U_a_b);
        	//potentialMaster.addPotential(U_a_b, new ISpecies[] {species,species} );

        // INTRAmolecular harmonic bending potential

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

        //      sim.mcMoveWiggle[0].setBondLength(1.54);
     //       sim.mcMoveWiggle[1].setBondLength(1.54);

        MCMoveClusterTorsionMulti[] torsionMoves = null;

        //  sim.mcMoveWiggle[0].setBondLength(1.54);
            //Ssim.mcMoveWiggle[1].setBondLength(1.54);

        //}

//         sim.integratorOS.setAdjustStepFreq(false);
//         sim.integratorOS.setStepFreq0(1);


        Box referenceBox = sim.box[0];
        Box targetBox = sim.box[1];

        sim.integratorOS.setNumSubSteps(1000); // Is this necessary?

        /*
        ****************************************************************************
        ****************************************************************************
        Directives for graphics

        true to run graphics (and not collect data)
        false to not run graphics (and collect data)
        ****************************************************************************
        ****************************************************************************
        */

        if (graphics) {

            referenceBox.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            targetBox.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(referenceBox).setShowBoundary(false);
            simGraphic.getDisplayBox(targetBox).setShowBoundary(false);

            simGraphic.getDisplayBox(referenceBox).setLabel("Reference-System Sampling");
            simGraphic.getDisplayBox(targetBox).setLabel("Target-System Sampling");

            // Create instances of ColorSchemeByType for reference and target simulations
            ColorSchemeByType colorScheme0 = (ColorSchemeByType) simGraphic.getDisplayBox(referenceBox).getColorScheme();
            ColorSchemeByType colorScheme1 = (ColorSchemeByType) simGraphic.getDisplayBox(targetBox).getColorScheme();

          /*  if (ethanol) {

            	SpeciesEthanol species = (SpeciesEthanol)sim.species;

            	// Create instances of the types of molecular sites
            	IAtomTypeLeaf type_O  = species.getOxygenType();
                IAtomTypeLeaf type_aC = species.getAlphaCarbonType();
                IAtomTypeLeaf type_C = species.getCarbonType();
                IAtomTypeLeaf type_aH = species.getAlphaHydrogenType();
                IAtomTypeLeaf type_H  = species.getHydrogenType();
                IAtomTypeLeaf type_X  = species.getXType();

                // Set color of each site type for each simulation

                colorScheme0.setColor(type_O, Color.RED);
                colorScheme0.setColor(type_aC, Color.GRAY);
                colorScheme0.setColor(type_C, Color.GRAY);
                colorScheme0.setColor(type_aH, Color.WHITE);
                colorScheme0.setColor(type_H, Color.WHITE);
                colorScheme0.setColor(type_X, Color.BLUE);

                colorScheme1.setColor(type_O, Color.RED);
                colorScheme1.setColor(type_aC, Color.GRAY);
                colorScheme1.setColor(type_C, Color.GRAY);
                colorScheme1.setColor(type_aH, Color.WHITE);
                colorScheme1.setColor(type_H, Color.WHITE);
                colorScheme1.setColor(type_X, Color.BLUE);

            }
            else {*/

            //SpeciesMethanol species = (SpeciesMethanol)sim.species;

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

            //}

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

        /*
        ****************************************************************************
        ****************************************************************************
        Other directives for simulation
        ****************************************************************************
        ****************************************************************************
        */

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+numMolecules+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);

        sim.setAccumulatorBlockSize((int)steps);

        System.out.println();
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        System.out.println();

        IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" blocks of 1000 attempted MC moves: ");
                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                double ratio = ratioAndError[0];
                double error = ratioAndError[1];
                System.out.println("Calculated B" + numMolecules + " = "+ratio*HSB[numMolecules]+" +/- "+error*HSB[numMolecules] + " Angstroms^3");

                DataGroup reference = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());

                System.out.println();
                System.out.println("Values calculated using the reference system's sampling");
                System.out.println("  average sign of reference system's integrand: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                        + "    stdev: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                        + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
                System.out.println("  average of overlap system's normalized integrand: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                        + "    stdev: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                        + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
                System.out.println("  ratio of these averages: " + ((DataDoubleArray) reference.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                        + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

                DataGroup targetData = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
                System.out.println();
                System.out.println("Values calculated using the target system's sampling");

                System.out.println("  average sign of target system's integrand: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                        + "    stdev: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                        + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
                System.out.println("  average of overlap system's normalized integrand: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                        + "    stdev: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                        + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
                System.out.println("  ratio of these averages: " + ((DataDoubleArray) targetData.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                        + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

                System.out.println();
                System.out.println("ratio calculated in target system divided by ratio calculated in reference system: "+ratio+", error: "+error);
                System.out.println("Calculated B" + numMolecules +  " = " +ratio*HSB[numMolecules]+" +/- "+error*HSB[numMolecules] + " Angstroms^3");
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval((int)(steps/10));
        sim.integratorOS.getEventManager().addListener(progressReportListener);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());


        DataGroup reference = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());

        System.out.println();
        System.out.println("Values calculated using the reference system's sampling");
        System.out.println("  average sign of reference system's integrand: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                + "    stdev: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
        System.out.println("  average of overlap system's normalized integrand: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + "    stdev: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
        System.out.println("  ratio of these averages: " + ((DataDoubleArray) reference.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + "    error: " + ((DataDoubleArray) reference.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        DataGroup targetData = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println();
        System.out.println("Values calculated using the target system's sampling");

        System.out.println("  average sign of target system's integrand: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                + "    stdev: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
        System.out.println("  average of overlap system's normalized integrand: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + "    stdev: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
        System.out.println("  ratio of these averages: " + ((DataDoubleArray) targetData.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + "    error: " + ((DataDoubleArray) targetData.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println();
        System.out.println("ratio calculated in target system divided by ratio calculated in reference system: "+ratio+", error: "+error);
        System.out.println("Calculated B" + numMolecules +  " = " +ratio*HSB[numMolecules]+" +/- "+error*HSB[numMolecules] + " Angstroms^3");
	}

    // to control whether or not graphics are used:
    
    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

        // number of molecules in simulation (e.g., 2 for B2 calculation)
    	public int numMolecules = 2;

        public double temperature = 400.0;   // Kelvin

        // number of overlap sampling steps
        // for each overlap sampling step, the simulation boxes are allotted
        // 1000 attempts for MC moves, total
        public long numSteps = 1000;


    }
    
}


