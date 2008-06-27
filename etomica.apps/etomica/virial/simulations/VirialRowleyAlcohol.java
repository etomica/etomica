package etomica.virial.simulations;

import java.awt.Color;

import etomica.api.IAction;
import etomica.api.IAtomTypeLeaf;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.rowley.EthanolPotentialHelper;
import etomica.models.rowley.MethanolPotentialHelper;
import etomica.models.rowley.SpeciesEthanol;
import etomica.models.rowley.SpeciesMethanol;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryEthanol;
import etomica.virial.SpeciesFactoryMethanol;
import etomica.virial.cluster.Standard;

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


    public static void main(String[] args) {
    	
        VirialParam params = new VirialParam();
        
        if (args.length > 0) {
            ReadParameters paramReader = new ReadParameters(args[0], params);
            paramReader.readParameters();
        }
        
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        boolean ethanol = params.ethanol;
        boolean pointCharges = params.pointCharges;
        
        // Diameter of hard spheres in reference system
        // Should be about the size of the molecules in the target system 
        double sigmaHSRef;
        if (ethanol) {
        	sigmaHSRef = 5;
        }
        else {
            sigmaHSRef = 3; 
        }
        
        final double[] HSB = new double[8];
        
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        
        System.out.println();
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+ " B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
		
        Space space = Space3D.getInstance();
           
             
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for overlap sampling
        ****************************************************************************
        ****************************************************************************
        */
        
        MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        
        // U_a_b is a pairwise potential (2 molecules, a and b, are involved).
        // The directives for calculation of U_a_b are provided later.
        PotentialGroup U_a_b = new PotentialGroup(2, space);
        
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
        System.out.println("B"+nPoints+" at "+temperature+"K");
         
        temperature = Kelvin.UNIT.toSim(temperature); // What are the simulation units for T?
         
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
         
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        
        final SimulationVirialOverlap sim;
        
        
        if(ethanol) {
        	sim = new SimulationVirialOverlap(space,new SpeciesFactoryEthanol(pointCharges),
        			temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesEthanol species = (SpeciesEthanol)sim.species;
        	EthanolPotentialHelper.initPotential(space, species, U_a_b, pointCharges);
        }
        else {
        	sim = new SimulationVirialOverlap(space,new SpeciesFactoryMethanol(pointCharges),
                    temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
        	SpeciesMethanol species = (SpeciesMethanol)sim.species;
        	MethanolPotentialHelper.initPotential(space, species, U_a_b, pointCharges);
        }
         
//         sim.integratorOS.setAdjustStepFreq(false);
//         sim.integratorOS.setStepFreq0(1);
        
         
        
        
        // Directives for calculating molecular potential, U_a_b
        // Cannot be called until species is instantiated.
        
        
        
       

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
             
        if (false) { 
            
            sim.box[0].getBoundary().setDimensions(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setDimensions(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            
            // Create instances of ColorSchemeByType for reference and target simulations
            ColorSchemeByType colorScheme0 = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme();
            ColorSchemeByType colorScheme1 = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            
            if (ethanol) {
            	
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
            else {
            	
            	SpeciesMethanol species = (SpeciesMethanol)sim.species;
            	
            	// Create instances of the types of molecular sites
            	IAtomTypeLeaf type_O  = species.getOxygenType();
                IAtomTypeLeaf type_aC = species.getAlphaCarbonType(); 
                IAtomTypeLeaf type_aH = species.getAlphaHydrogenType();
                IAtomTypeLeaf type_H  = species.getHydrogenType();
                IAtomTypeLeaf type_X  = species.getXType();
                
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
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
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
        
        IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double ratio = sim.dsvo.getDataAsScalar();
                double error = sim.dsvo.getError();
                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
            }
        };
        sim.integratorOS.addIntervalAction(progressReport);
        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        
        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);

        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        
        System.out.println();
        System.out.println("reference ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("reference average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("reference overlap function average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println();
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("target average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("target overlap function average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
    	
    	// number of molecules in simulation (e.g., 2 for B2 calculation)
    	public int nPoints = 2;
        
        public double temperature = 500.0;   // Kelvin
        
        public long numSteps = 10000;
        
        // ethanol = false: methanol
        // ethanol = true: ethanol
        protected boolean ethanol = true;
        
        // pointCharges = false: Rowley et al (2006) model without point charges
        // pointCharges = true: Rowley et al (2006) model with point charges
        protected boolean pointCharges = true;
    }
    
}

