/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.action.IntegratorDimerApproach;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataLogger.DataWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.*;

import java.awt.*;

/*
 * Uses the site-site models of Rowley et al (2006) to reproduce potential-energy plots for dimers 
 * of ethanol or methanol molecules along particular approach routes (please see the paper for details).
 * 
 * Set the last three fields at the end of the file to specify the species, the model (with or without
 * point charges), and the approach route.  
 * 
 * Please note that about half of the approach routes do not generate the correct initial configuration, and 
 * the implementation of the potentials with point charges has yet to be validated.  In other words, this
 * class is still junk.
 * 
 * K.R. Schadel 2008 
 */


public class DimerApproach extends Simulation {
	
    // True to use model with point charges, false to use model without point charges
    public final static boolean pointCharges = true;
	public final static long serialVersionUID = 1L;
	public static SpeciesEthanol speciesEthanol;
	public static SpeciesMethanol speciesMethanol;
	public static MeterPotentialEnergy meterPE;
    // True to consider ethanol, false to consider methanol
    static boolean ethanol = false;
    // ID of approach route (see Rowley et al (2006) for table)
    static int route = 18;
    static IAtom atom_O_A;
    static IAtom atom_aC_A;
    static IAtom atom_aH_A;
    static IAtom atom_H1_A;
    static IAtom atom_O_B;
    static IAtom atom_aC_B;
    static IAtom atom_aH_B;
    public ISpecies species;
    public Box box;
    public PotentialMaster potentialMaster;
    public IntegratorDimerApproach dimerApproach;
    
    

	public DimerApproach() {


        super(Space3D.getInstance());

        // *************
        // The Species
        // *************
        if (ethanol) {
            species = new SpeciesEthanol(space, pointCharges);
            speciesEthanol = (SpeciesEthanol) species;

        } else {
            species = new SpeciesMethanol(space, pointCharges);
            speciesMethanol = (SpeciesMethanol) species;
        }
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));

        // ***************
        // The Potential
        // ***************

        PotentialGroup U_a_b = new PotentialGroup(2);
        if (ethanol) {
            EthanolPotentialHelper.initPotential(space, speciesEthanol, U_a_b, pointCharges);

        } else {
            double sigmaOC = 0.00001;
            double sigmaOH = 0.05;
            MethanolPotentialHelper.initPotential(space, speciesMethanol, U_a_b, pointCharges, sigmaOC, sigmaOH);
        }
        addSpecies(species);
        box.setNMolecules(species, 2); // 2 molecules in box...
        U_a_b.setBox(box);
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(U_a_b, new ISpecies[]{species, species});

        // *********************
        // The Integrator
        // *********************

        dimerApproach = new IntegratorDimerApproach(potentialMaster, box);

        // Methods in dimerApproach that must be called
        dimerApproach.setMolecules();
        dimerApproach.setImportantAtoms();
        dimerApproach.setRoute(route);

        double[][] params;
        if (ethanol) {
            params = EthanolRouteParams.setEthanolParams(route);
        } else {
            params = MethanolRouteParams.setMethanolParams(route);
        }

        dimerApproach.setRouteParams(params);

        // The following is required for dataDistances:
        atom_O_A = dimerApproach.getAtom_O_A();
        atom_aC_A = dimerApproach.getAtom_aC_A();
        atom_aH_A = dimerApproach.getAtom_aH_A();
        atom_H1_A = dimerApproach.getAtom_H1_A();

        atom_O_B = dimerApproach.getAtom_O_B();
        atom_aC_B = dimerApproach.getAtom_aC_B();
        atom_aH_B = dimerApproach.getAtom_aH_B();


        // This may be called here or in the main method
        // Must be called after above set methods
        dimerApproach.initializeCoordinates();

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(dimerApproach);
        activityIntegrate.setMaxSteps(75);
        activityIntegrate.setSleepPeriod(100);

        getController().addAction(activityIntegrate);

    }

	public static void main(String[] string)  {
		
		DimerApproach sim = new DimerApproach();
		
		/*
	     ****************************************************************************
	     ****************************************************************************
         Directives for calculating and storing potential energy
         ****************************************************************************
	     ****************************************************************************
	     */

		meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
		
		DataLogger dataLoggerPE = new DataLogger();
		
		DataFork dataForkPE = new DataFork();

		DataPump dataPumpPE = new DataPump(meterPE, dataForkPE);
		
		dataForkPE.addDataSink(dataLoggerPE);
		
		sim.dimerApproach.getEventManager().addListener(new IntegratorListenerAction(dataPumpPE)); // measure data at each step
		
		dataLoggerPE.setFileName("Potential energy");
		DataWriter dataWriterR = new DataTableWriter();
        dataLoggerPE.setDataSink(dataWriterR);

		/*
	     ****************************************************************************
	     ****************************************************************************
         Directives for calculating and storing distances  between sites
         ****************************************************************************
	     ****************************************************************************
	     */
	   
		String label1 = "Distance between alpha carbons (Angstroms) ";
		String label2 = "Distance between oxygen (monomer B) and alpha hydrogen (monomer A) (Angstroms) ";
		String label3 = "Distance between oxygen (monomer A) and alpha hydrogen (monomer B) (Angstroms) (Angstroms) ";
		String label4 = "Distance between oxygens (Angstroms) ";
		String labelGrr = "Distance between oxygen (monomer B) and hydrogen1 (monomer A) (Angstroms) ";
		
        DataSourceAtomDistance  dataDistance1 = new DataSourceAtomDistance(sim.space);
        DataSourceAtomDistance  dataDistance2 = new DataSourceAtomDistance(sim.space);
        DataSourceAtomDistance  dataDistance3 = new DataSourceAtomDistance(sim.space);
        DataSourceAtomDistance  dataDistance4 = new DataSourceAtomDistance(sim.space);
        DataSourceAtomDistance  dataDistanceGrr = new DataSourceAtomDistance(sim.space);
        
        
		dataDistance1.setAtoms(atom_aC_A, atom_aC_B);
		dataDistance2.setAtoms(atom_O_B,  atom_aH_A);
		dataDistance3.setAtoms(atom_O_A,  atom_aH_B);
		dataDistance4.setAtoms(atom_O_A,  atom_O_B);
		dataDistanceGrr.setAtoms(atom_O_B,  atom_H1_A);
		
		
	/*
		DataLogger dataLoggerDistance1 = new DataLogger();
		DataLogger dataLoggerDistance2 = new DataLogger();
		DataLogger dataLoggerDistance3 = new DataLogger();

		DataPump dataPumpDistance1 = new DataPump(dataDistance1, dataLoggerDistance1);
		DataPump dataPumpDistance2 = new DataPump(dataDistance2, dataLoggerDistance2);
		DataPump dataPumpDistance3 = new DataPump(dataDistance3, dataLoggerDistance3);
		
		sim.dimerApproach.addIntervalAction(dataPumpDistance1); // measure data at each step
		sim.dimerApproach.addIntervalAction(dataPumpDistance2); 
		sim.dimerApproach.addIntervalAction(dataPumpDistance3); 

	*/
		/*dataLoggerDistance1.setFileName("Distance between alpha carbons");
		dataLoggerDistance2.setFileName("Distance between oxygen (monomer B) and alpha hydrogen (monomer A)");
		dataLoggerDistance3.setFileName("Distance between oxygen (monomer A) and alpha hydrogen (monomer B)");
		*/
		
		/*
	     ****************************************************************************
	     ****************************************************************************
         Directives for graphics
	        true to run graphics 
            false to not run graphics 
         ****************************************************************************
	     ****************************************************************************
	     */

        if (true) { 
            
            sim.box.getBoundary().setBoxSize(sim.space.makeVector(new double[]{40,40,40}));
            
            // *********************
            // The Title
            // *********************
            
            String string1;
            String string2;
            
            if (ethanol) {
            	string1 = "Ethanol ";
            } else {
            	string1 = "Methanol ";
            }
            
            if (pointCharges) {
            	string2 = "with point charges";
            } else {
            	string2 = "without point charges";
            }
            
            String appName = string1 + string2 +  ": Route " + route;
            
            // ****************************************
            // Things that matter more than the title
            // ****************************************
            
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName);
            // The default Paint Interval is too infrequent
            simGraphic.setPaintInterval(sim.box, 1);
            //simGraphic.getDisplayBox(sim.box).setShowBoundary(false);
            //simGraphic.getDisplayBox(sim.box).setLabel("ARGGGG");
  
            
            // Create instances of ColorSchemeByType for reference and target simulations
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme();
	            
        	
            if (ethanol) {
            	
            	// Create instances of the types of molecular sites

                AtomType type_O = speciesEthanol.getOxygenType();
                AtomType type_aC = speciesEthanol.getAlphaCarbonType();
                AtomType type_C = speciesEthanol.getCarbonType();
                AtomType type_aH = speciesEthanol.getAlphaHydrogenType();
                AtomType type_H = speciesEthanol.getHydrogenType();
                AtomType type_X = speciesEthanol.getXType();
                
                // Set color of each site type for each simulation
                
                colorScheme.setColor(type_O, Color.RED);
                colorScheme.setColor(type_aC, Color.GRAY);
                colorScheme.setColor(type_C, Color.GRAY);
                colorScheme.setColor(type_aH, Color.WHITE);
                colorScheme.setColor(type_H, Color.WHITE);
                colorScheme.setColor(type_X, Color.BLUE);
	        	
            } else {
            	
            	// Create instances of the types of molecular sites

                AtomType type_O = speciesMethanol.getOxygenType();
                AtomType type_aC = speciesMethanol.getAlphaCarbonType();
                AtomType type_aH = speciesMethanol.getAlphaHydrogenType();
                AtomType type_H = speciesMethanol.getHydrogenType();
                AtomType type_X = speciesMethanol.getXType();
                
                // Set color of each site type for each simulation
                
                colorScheme.setColor(type_O, Color.RED);
                colorScheme.setColor(type_aC, Color.GRAY);
                colorScheme.setColor(type_aH, Color.WHITE);
                colorScheme.setColor(type_H, Color.WHITE);
                colorScheme.setColor(type_X, Color.BLUE);
            	
            }
           
            /*
    	     ****************************************************************************
    	     ****************************************************************************
             Plotting potential energy vs. distance between alpha carbons
             ****************************************************************************
    	     ****************************************************************************
    	     */
            
        
            AccumulatorHistory energyHistory = new AccumulatorHistory();
            
            // To use distance between alpha carbons as independent variable for plot, rather than the step number
            energyHistory.setTimeDataSource(dataDistance1);
           
            dataForkPE.addDataSink(energyHistory);
            
            DisplayPlot ePlot = new DisplayPlot();
            
            ePlot.setLabel("Potential Energy");
            
            // Create conversion factor to change epsilon values from simulation units to kcal/mol
    		// NOTE: This only works for the graphics, not the dataLogger
    		Unit eUnit = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
            ePlot.setUnit(eUnit);
            
            energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
            
    		//ePlot.setDoLegend(true);
    		
    		ePlot.getPlot().setXRange(-2.0, 50.0);
    		
    		if (ethanol) {
    			if (route == 15 || route == 16) {
    				ePlot.getPlot().setYRange(-6.0, 2.0);
    			} else {
    				ePlot.getPlot().setYRange(-2.5, 2.0);
    			}
    		} else {
    			if (route >= 11 && route <= 17) {
    				ePlot.getPlot().setYRange(-6, 2.0);
    			} else {
    				ePlot.getPlot().setYRange(-2.0, 2.0);
    			}
    		}
    		
    		ePlot.getPlot().setXLabel("Distance between alpha carbons (Angstroms)");
    		ePlot.getPlot().setYLabel("Potential energy (kcal/mol)");
    		
    		simGraphic.add(ePlot);
    		
    		/*
    	     ****************************************************************************
    	     ****************************************************************************
             Plotting potential energy vs. distance between 
             oxygen(monomer B) and hydrogen(2) (monomer A)
             ****************************************************************************
    	     ****************************************************************************
    	     */
            
        
           /* AccumulatorHistory energyHistoryGrr = new AccumulatorHistory();
            
            // To use distance between alpha carbons as independent variable for plot, rather than the step number
            energyHistoryGrr.setTimeDataSource(dataDistanceGrr);
           
            dataForkPE.addDataSink(energyHistoryGrr);
            
            DisplayPlot ePlotGrr = new DisplayPlot();
            
            ePlotGrr.setLabel("Potential Energy vs. rOH");
            
            energyHistoryGrr.setDataSink(ePlotGrr.getDataSet().makeDataSink());
            
    		ePlotGrr.setDoLegend(true);
    		
    		simGraphic.add(ePlotGrr);*/
            
            simGraphic.makeAndDisplayFrame(appName);
            
            return;
            
        } else {
        	// this method is called when graphics are used and the start button is pressed 
        	// must be included here:
    		sim.getController().actionPerformed();
        }
	}
	
	
	
}
