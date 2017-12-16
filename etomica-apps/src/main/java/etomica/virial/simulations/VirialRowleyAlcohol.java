/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.DataLogger.DataWriter;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.models.rowley.*;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Mayer-sampling MC simulation for methanol or ethanol using the potentials published in:
 *   R.L. Rowley, C.M. Tracy, & T.A. Pakkanen. (2006)  
 *   "Potential energy surfaces for small alcohol dimers I: Methanol and ethanol."  
 *   J. Chem. Phys.  125. 154302. 
 *   
 * Set species, model (with or without point charges), virial-coefficent type (B2, B3, etc.),
 * and temperature in VirialParam() method at bottom of file. 
 *   
 * Class adapted from VirialAlkane by K.R. Schadel, May 2008
 *  
 */


// ****************************************************************************
// About the potential: 
// ****************************************************************************
  /*
  The Morse potential (between sites on different molecules) is of the form:
  
     u(i,j) = -epsilon(i,j) * ( 1 - { 1 - exp( -A(i,j) ( r(i,j) - re(i,j) ) }^2 ) 

  The molecular potential (between molecules a and b) is of the form:
  
     U(a,b) = sum over i (sites of a) { sum over j (sites of b) { u(i,j) }}

  Because bond lengths and angles are fixed, we do not need to compute INTRAmolecular site-site interactions.
  
  The model with point charges is modified to include a Coulombic interaction.
  */
// ****************************************************************************


public class VirialRowleyAlcohol {

    public static MeterPotentialEnergy meterPE;
    // ethanol = false: methanol
    // ethanol = true: ethanol
    protected static boolean ethanol = false;
    // pointCharges = false: Rowley et al (2006) model without point charges
    // pointCharges = true: Rowley et al (2006) model with point charges
    protected static boolean pointCharges = true;
    protected static boolean printapalooza = true;

    // to control whether or not graphics are used:
    protected static boolean graphics = false;
    protected static boolean plots = false;
    protected static double sigmaOC = 0.00001;
    protected static double sigmaOH = 0.05;
    
    public static void main(String[] args) {

        VirialParam params = new VirialParam();

        /*System.out.println(""+ (Double.POSITIVE_INFINITY-Double.POSITIVE_INFINITY));
        System.out.println(""+ Math.sqrt(Double.POSITIVE_INFINITY));
        System.out.println("" + (Double.POSITIVE_INFINITY)*(Double.POSITIVE_INFINITY));
        System.exit(1);*/

        // enables one to overwrite parameters values in VirialParam() and use those provided in string instead

     // enables one to overwrite parameters values in VirialParam() and use those provided in string instead
        final int numMolecules;
        double temperature;
        long steps;
        if (args.length == 0) {
        	numMolecules = params.numMolecules;
            temperature = params.temperature;

            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            steps = params.numSteps;
        } else if (args.length == 4 || args.length == 3) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	numMolecules = Integer.parseInt(args[0]);
        	temperature = Integer.parseInt(args[1]);
            steps = Integer.parseInt(args[2]);

        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed to VirialRowleyAlcohol.");
        }


        // Diameter of hard spheres in reference system
        // Should be about the size of the molecules in the target system
        double sigmaHSRef;
        if (ethanol) {
        	sigmaHSRef = 5;
        }
        else {
            sigmaHSRef = 6;
        }

        final double[] HSB = new double[8];

        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);

        if (printapalooza) {

	        System.out.println();
	        System.out.println("sigmaHSRef: "+sigmaHSRef + " Angstroms");
	        System.out.println("B2HS: "+HSB[2] + " Angstroms^3");
	        System.out.println("B3HS: "+HSB[3] + " Angstroms^6"); // + " = " +(HSB[3]/(HSB[2]*HSB[2]))+ " B2HS^2");
	        System.out.println("B4HS: "+HSB[4] + " Angstroms^9"); // + " = " +(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
	        System.out.println("B5HS: "+HSB[5] + " Angstroms^12"); // What's with this?: + " = 0.110252 B2HS^4");
	        System.out.println("B6HS: "+HSB[6] + " Angstroms^15"); //  + " = 0.03881 B2HS^5");
	        System.out.println("B7HS: "+HSB[7] + " Angstroms^18"); // + " = 0.013046 B2HS^6");
        }

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
        System.out.print("Rowley et al model for");
        if (ethanol) {
        	System.out.print(" ethanol");
        }
        else {
        	System.out.print(" methanol");
        }
        if (pointCharges) {
        	System.out.print(" with point charges");
        }
        else {
        	System.out.print(" without point charges");
        }
        System.out.println();
        System.out.println("B"+numMolecules+" at "+temperature+"K");

        temperature = Kelvin.UNIT.toSim(temperature); // What are the simulation units for T?

        ClusterAbstract targetCluster = Standard.virialCluster(numMolecules, fTarget, numMolecules>3, eTarget, true);

        if (pointCharges) {
            ((ClusterSum)targetCluster).setCaching(false);
        	targetCluster = new ClusterCoupledFlipped(targetCluster, space);
        }
        targetCluster.setTemperature(temperature);

        ClusterAbstract refCluster = Standard.virialCluster(numMolecules, fRef, numMolecules>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println(steps*1000+" total attempted MC moves ("+steps+" blocks of 1000)");

        final SimulationVirialOverlap sim;
        PotentialMaster potentialMaster = new PotentialMaster();

        if(ethanol) {
        	sim = new SimulationVirialOverlap (space,new SpeciesFactoryEthanol(pointCharges),
        			temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesEthanol species = (SpeciesEthanol)sim.getSpecies(0);
        	EthanolPotentialHelper.initPotential(space, species, U_a_b, pointCharges);
        	//potentialMaster.addPotential(U_a_b, new ISpecies[] {species,species} );
        }
        else {
        	sim = new SimulationVirialOverlap (space,new SpeciesFactoryMethanol(pointCharges),
                    temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesMethanol species = (SpeciesMethanol)sim.getSpecies(0);

        	MethanolPotentialHelper.initPotential(space, species, U_a_b, pointCharges, sigmaOC, sigmaOH);
        	//potentialMaster.addPotential(U_a_b, new ISpecies[] {species,species} );
        }

//         sim.integratorOS.setAdjustStepFreq(false);
//         sim.integratorOS.setStepFreq0(1);


        BoxCluster referenceBox = sim.box[0];
        BoxCluster targetBox = sim.box[1];


        sim.integratorOS.setNumSubSteps(1000); // Is this necessary?

        /*
         ****************************************************************************
         ****************************************************************************
         Set the seed
         ****************************************************************************
         ****************************************************************************
         */

        if (args.length == 4 ) {
        	long seed = Integer.parseInt(args[3]);
        	//((RandomNumberGenerator) sim.random).getWrappedRandom().setSeed(seed);
        	System.out.println();
        	System.out.println("Trying seed = "+ seed);
        	System.out.println();
        }

        /*
        ****************************************************************************
        ****************************************************************************
        Directives for graphics

        true to run graphics (and not collect data)
        false to not run graphics (and collect data)
        ****************************************************************************
        ****************************************************************************
        */

        if (plots) {


        	/* *************************************************************
             * Separation Distance Histogram
             * *************************************************************
             */

            String label1 = "Distance between alpha carbons (Angstroms) ";

            DataSourceAtomDistance  dataDistance1 = new DataSourceAtomDistance(space);
    		DataFork dataForkDistance = new DataFork();
    		DataPump dataPumpDistance = new DataPump(dataDistance1, dataForkDistance);

            sim.integrators[1].getEventManager().addListener(new IntegratorListenerAction(dataPumpDistance)); // measure data at each step

            IMoleculeList moleculeList = targetBox.getMoleculeList();
    		IMolecule monomerA = moleculeList.getMolecule(0);
    		IMolecule monomerB = moleculeList.getMolecule(1);

            IAtomList atomSetA = monomerA.getChildList();
    		IAtomList atomSetB = monomerB.getChildList();

            IAtom atom_aC_A = atomSetA.getAtom(1);
    		IAtom atom_aC_B = atomSetB.getAtom(1);

            dataDistance1.setAtoms(atom_aC_A, atom_aC_B);

            AccumulatorHistogram accumulatorR = new AccumulatorHistogram();
            dataForkDistance.addDataSink(accumulatorR);

            DataLogger dataLoggerR = new DataLogger();
            dataLoggerR.setFileName("Separation Distance");
            DataWriter dataWriterR = new DataTableWriter();
            dataLoggerR.setDataSink(dataWriterR);
            dataLoggerR.setAppending(false);
            dataLoggerR.setCloseFileEachTime(true);
            accumulatorR.setDataSink(dataLoggerR);
            accumulatorR.setPushInterval(100);

    		/* *************************************************************
             * Sampling Weight vs. Separation Distance
             * *************************************************************
             */

            MeterSamplingWeight meterPi = new MeterSamplingWeight(dataDistance1);
    		meterPi.setBox(targetBox);

            DataPump dataPumpPi = new DataPump(meterPi, null);

            DataLogger dataLoggerPi = new DataLogger();


            DataFork dataForkPi = new DataFork();

    		DataPump dataPumpPi2 = new DataPump(meterPi, dataForkPi);

            sim.integrators[1].getEventManager().addListener(new IntegratorListenerAction(dataPumpPi));
    		sim.integrators[1].getEventManager().addListener(new IntegratorListenerAction(dataPumpPi2));

            //dataForkPi.addDataSink(dataLoggerPi);

            //dataForkDistance.addDataSink(dataLoggerR);

            dataLoggerPi.setFileName("Sampling Weight");

            //

            HistogramNotSoSimple histogramPi = new HistogramNotSoSimple(12000, new DoubleRange(0,1200));
            AccumulatorHistogram accumulatorPi = new AccumulatorHistogram(histogramPi);
            dataPumpPi.setDataSink(accumulatorPi);

            DataWriter dataWriterPi = new DataTableWriter();
            dataLoggerPi.setDataSink(dataWriterPi);
            dataLoggerPi.setAppending(false);
            dataLoggerPi.setCloseFileEachTime(true);

            accumulatorPi.addDataSink(dataLoggerPi);

            accumulatorPi.setPushInterval(100);


        }

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

            if (ethanol) {

                SpeciesEthanol species = (SpeciesEthanol)sim.getSpecies(0);

                // Create instances of the types of molecular sites
                AtomType type_O = species.getOxygenType();
                AtomType type_aC = species.getAlphaCarbonType();
                AtomType type_C = species.getCarbonType();
                AtomType type_aH = species.getAlphaHydrogenType();
                AtomType type_H = species.getHydrogenType();
                AtomType type_X = species.getXType();

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
            else {

                SpeciesMethanol species = (SpeciesMethanol)sim.getSpecies(0);

                // Create instances of the types of molecular sites
                AtomType type_O = species.getOxygenType();
                AtomType type_aC = species.getAlphaCarbonType();
                AtomType type_aH = species.getAlphaHydrogenType();
                AtomType type_H = species.getHydrogenType();
                AtomType type_X = species.getXType();

                // Set color of each site type for each simulation

                colorScheme0.setColor(type_O, Color.RED);
                colorScheme0.setColor(type_aC, Color.GRAY);
                colorScheme0.setColor(type_aH, Color.WHITE);
                colorScheme0.setColor(type_H, Color.WHITE);
                colorScheme0.setColor(type_X, Color.BLUE);

                colorScheme1.setColor(type_O, Color.RED);
                colorScheme1.setColor(type_aC, Color.GRAY);
                colorScheme1.setColor(type_aH, Color.WHITE);
                colorScheme1.setColor(type_H, Color.WHITE);
                colorScheme1.setColor(type_X, Color.BLUE);
            }
        /*
            if (plots) {


            	 *************************************************************
                 * Separation Distance Histogram
                 * *************************************************************





                DisplayPlot rPlot = new DisplayPlot();

                rPlot.setLabel("Separation Distance");
                rPlot.setDoLegend(false);

                accumulatorR.setDataSink(rPlot.getDataSet().makeDataSink());
                accumulatorR.setPushInterval(1000);

                simGraphic.add(rPlot);

	    		 *************************************************************
                 * Potential Energy Plot
                 * *************************************************************


	    		meterPE = new MeterPotentialEnergy(potentialMaster);

        		meterPE.setBox(targetBox);

        		DataLogger dataLoggerPE = new DataLogger();

        		DataFork dataForkPE = new DataFork();

        		DataPump dataPumpPE = new DataPump(meterPE, dataForkPE);

        		DataPump dataPumpPE2 = new DataPump(meterPE, null);

        		dataForkPE.addDataSink(dataLoggerPE);

        		dataLoggerPE.setFileName("Potential energy");

        		sim.integrators[1].addIntervalAction(dataPumpPE);

                HistogramNotSoSimple histogramPE = new HistogramNotSoSimple(1000, new DoubleRange(0,1200));
                AccumulatorHistogram accumulatorPE = new AccumulatorHistogram(histogramPE);
                dataForkPE.setDataSink(accumulatorPE);

                DisplayPlot ePlot = new DisplayPlot();

                ePlot.setLabel("Potential Energy");
                ePlot.setDoLegend(false);
                ePlot.getPlot().setXRange(3.0, 1200);

                accumulatorPE.setDataSink(ePlot.getDataSet().makeDataSink());
                accumulatorPE.setPushInterval(1000);

                simGraphic.add(ePlot);


                 *************************************************************
                 * Sampling Weight vs. Separation Distance
                 * *************************************************************



                HistogramNotSoSimple histogramPi = new HistogramNotSoSimple(1000, new DoubleRange(0,1200));
                AccumulatorHistogram accumulatorPi = new AccumulatorHistogram(histogramPi);
                dataPumpPi.setDataSink(accumulatorPi);

                DisplayPlot piPlot = new DisplayPlot();

                piPlot.setLabel("Sampling Weight");
                piPlot.setDoLegend(false);

                accumulatorPi.setDataSink(piPlot.getDataSet().makeDataSink());
                accumulatorPi.setPushInterval(1000);

                simGraphic.add(piPlot);



                /* *************************************************************
                 * Paint Interval
                 * *************************************************************




            }

            */

            simGraphic.setPaintInterval(targetBox, 1);

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

        /*System.out.println();
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        System.out.println();*/

        IAction progressReport = new IAction() {
            public void actionPerformed() {

                if (printapalooza) {
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
                    System.out.println("Calculated B" + numMolecules + " = " + ratio * HSB[numMolecules] + " +/- " + error * HSB[numMolecules] + " Angstroms^3");
                    System.out.println();
            	}
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval((int)(steps/10));
        sim.integratorOS.getEventManager().addListener(progressReportListener);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        if (printapalooza) {
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
        }

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];

        if (printapalooza) {
	        System.out.println();
	        System.out.println("ratio calculated in target system divided by ratio calculated in reference system: "+ratio+", error: "+error);
	        System.out.println("Calculated B" + numMolecules +  " = " +ratio*HSB[numMolecules]+" +/- "+error*HSB[numMolecules] + " Angstroms^3");
        } else {
        	double coeff = ratio*HSB[numMolecules];
	        if (coeff > 0) {
	        	long seed = Integer.parseInt(args[3]);
	        	System.out.println("*******************************************************");
	        	System.out.println("B"+ numMolecules + " = " + coeff + "for seed = " + seed);
	        	System.out.println("*******************************************************");
	        	System.out.println();
	        	System.out.println();
	        	System.out.println();
            }
        }
	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

        // number of molecules in simulation (e.g., 2 for B2 calculation)
    	public int numMolecules = 2;

        public double temperature = 700.0;   // Kelvin

        // number of overlap sampling steps
        // for each overlap sampling step, the simulation boxes are allotted
        // 1000 attempts for MC moves, total
        public long numSteps = 10000;

        public long seed = 1;


    }
    
}

