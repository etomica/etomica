package etomica.virial.simulations;

import etomica.api.IAction;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.atom.iterator.ApiBuilder;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.models.rowley.P2RepRowley;
import etomica.models.rowley.SpeciesMethanol;
import etomica.potential.P2ModifiedMorse;
import etomica.potential.P2Morse;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Calorie;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.UnitRatio;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryMethanol;
import etomica.virial.cluster.Standard;

/**
 * Mayer sampling simulation for methanol using the multi-site modified Morse potential published in:
 *   R.L. Rowley, C.M. Tracy, & T.A. Pakkanen. (2006)  
 *   "Potential energy surfaces for small alcohol dimers I: Methanol and ethanol."  
 *   J. Chem. Phys.  125. 154302. 
 *   
 *  Class adapted from VirialAlkane by K.R. Schadel, May 2008
 *  
 */


// ****************************************************************************
// About the potential: 
// ****************************************************************************
  /*
  The modified Morse potential (between sites on different molecules) is of the form:
  
     u(i,j) = -epsilon(i,j) * ( 1 - { 1 - exp( -A(i,j) ( r(i,j) - re(i,j) ) }^2 ) + z(i)*z(j)*(electron_charge)^2/(4*pi*free_space_permittivity*r(i,j))

  The molecular potential (between molecules a and b) is of the form:
  
     U(a,b) = sum over i (sites of a) { sum over j (sites of b) { u(i,j) }}

  Because bond lengths and angles are fixed, we do not need to compute INTRAmolecular site-site interactions.
  */
// ****************************************************************************


public class VirialRowleyMethanolWPC {


    public static void main(String[] args) {
    	
        VirialParam params = new VirialParam();
        
        if (args.length > 0) {
            ReadParameters paramReader = new ReadParameters(args[0], params);
            paramReader.readParameters();
        }
        final int nPoints = params.nPoints;

        double temperature = params.temperature;
        long steps = params.numSteps;
        
        double sigmaHSRef = 3; //Initial guess: anything larger than molecular diameters in target system
        
        final double[] HSB = new double[8];
        
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        
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
        Parameters for site-site interactions between methanol molecules
        - Rowley model WITH point charges (modified Morse potential):
          
        The original units are:
          
        	kcal/mol for epsilon(i,j) - not in simulation units
        	1/Angstroms for A(i,j) - in simulation units
        	Angstroms for r and re - in simulation units
        	
       	Note: P2ModifiedMorse takes alpha, a dimensionalness parameter, rather than A
       	
        ****************************************************************************
        ****************************************************************************
        */
        
        // Create conversion factor to change epsilon values from kcal/mol to simulation units
        UnitRatio eunit = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
        
        // oxygen and oxygen (O-O)
        double epsilon_O_O = eunit.toSim(0.28128);
        double A_O_O = 2.01169;
        double re_O_O = 3.14862;
        
        // oxygen and the "alpha" carbon (O-aC)
        double epsilon_O_aC = eunit.toSim(0.00049);
        double A_O_aC = 2.87970;
        double re_O_aC = 3.98955;
        
        // oxygen and the "alpha" hydrogen (O-aH)
        double epsilon_O_aH = eunit.toSim(0.00953);
        double A_O_aH = 2.39278;
        double re_O_aH = 0.96556;
        
        // oxygen and hydrogen (O-H)
        double epsilon_O_H = eunit.toSim(0.00000479);
        double A_O_H = 1.44640;
        double re_O_H = 6.75499;
        
        // "alpha" carbon and "alpha" carbon(aC-aC)
        double epsilon_aC_aC = eunit.toSim(0.53431);
        double A_aC_aC = 5.42872;
        double re_aC_aC = 0.54673;
        
        // "alpha" carbon and "alpha" hydrogen (aC-aH)
        double epsilon_aC_aH = eunit.toSim(1.20290);
        double A_aC_aH = 1.48559;
        double re_aC_aH = 2.41005;
        
        // "alpha" carbon and hydrogen (aC-H)
        double epsilon_aC_H = eunit.toSim(0.30649);
        double A_aC_H = 1.80093;
        double re_aC_H = 2.67419;
        
        // "alpha" hydrogen and "alpha" hydrogen (aH-aH)
        double epsilon_aH_aH = eunit.toSim(0.00303);
        double A_aH_aH = 1.59429;
        double re_aH_aH = 3.81180;
        
        // "alpha" hydrogen and hydrogen (aH-H)
        double epsilon_aH_H = eunit.toSim(0.0000908);
        double A_aH_H = 0.18342;
        double re_aH_H = 1.38361;
        
        // fudge site and fudge site (X-X)
        // Potential simplified to purely repulsive form: uXX = BXX * exp(-CXX*rXX)
        double BXX = eunit.toSim(0.28457);
        double CXX = 6.79460;
        double rOX = 0.09310; // just used to draw X site
        
        // "alpha" hydrogen and fudge site (aH-X)
        double epsilon_aH_X = eunit.toSim(0.03155);
        double A_aH_X = 1.81525;
        double re_aH_X = 2.98180;
        
        // hydrogen and hydrogen (H-H)
        double epsilon_H_H = eunit.toSim(0.01048);
        double A_H_H = 1.26072;
        double re_H_H = 3.97536;
        
        // Point charges
        double z_O = -0.6423;
        double z_aC = 0.2551;
        double z_aH = 0.3873;
        
        /* 
         * The point charges for the non-alpha hydrogens and the satellite site are zero: 
         * One can use the unmodified Morse potential for interactions involving these sites.
         */
        
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for calculation of site-site interaction energies:
        ****************************************************************************
        ****************************************************************************
        */
        
        P2ModifiedMorse u_O_O   = new P2ModifiedMorse(space, epsilon_O_O,   re_O_O,   A_O_O,   z_O,  z_O   );
        P2ModifiedMorse u_O_aC  = new P2ModifiedMorse(space, epsilon_O_aC,  re_O_aC,  A_O_aC,  z_O,  z_aC  );
        P2ModifiedMorse u_O_aH  = new P2ModifiedMorse(space, epsilon_O_aH,  re_O_aH,  A_O_aH,  z_O,  z_aH  );
        P2Morse         u_O_H   = new         P2Morse(space, epsilon_O_H,   re_O_H,   A_O_H   );
        
        P2ModifiedMorse u_aC_aC = new P2ModifiedMorse(space, epsilon_aC_aC, re_aC_aC, A_aC_aC, z_aC, z_aC );
        P2ModifiedMorse u_aC_aH = new P2ModifiedMorse(space, epsilon_aC_aH, re_aC_aH, A_aC_aH, z_aC, z_aH );
        P2Morse         u_aC_H  = new         P2Morse(space, epsilon_aC_H,  re_aC_H,  A_aC_H  );
        
        P2ModifiedMorse u_aH_aH = new P2ModifiedMorse(space, epsilon_aH_aH, re_aH_aH, A_aH_aH, z_aH, z_aH );
        P2Morse         u_aH_H  = new         P2Morse(space, epsilon_aH_H,  re_aH_H,  A_aH_H  );
        P2Morse         u_aH_X  = new         P2Morse(space, epsilon_aH_X,  re_aH_X,  A_aH_X  );
        
        P2Morse         u_H_H   = new         P2Morse(space, epsilon_H_H,   re_H_H,   A_H_H   );

        P2RepRowley u_X_X = new P2RepRowley (space, BXX, CXX);

        
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for overlap sampling
        ****************************************************************************
        ****************************************************************************
        */
         
        // U_a_b is a pairwise potential (2 molecules, a and b, are involved).
        PotentialGroup U_a_b = new PotentialGroup(2, space); 
        
        MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
         
        MayerGeneral fTarget = new MayerGeneral(U_a_b);
        MayerEGeneral eTarget = new MayerEGeneral(U_a_b);
         
        System.out.println("B"+nPoints+" at "+temperature+"K");
         
        temperature = Kelvin.UNIT.toSim(temperature); // What are the simulation units for T?
         
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
         
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        
        boolean pointCharges = true;
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactoryMethanol(pointCharges),
                           temperature,refCluster,targetCluster); //use first constructor; no need for intramolecular movement MC trial
         
//         sim.integratorOS.setAdjustStepFreq(false);
//         sim.integratorOS.setStepFreq0(1);
        
        /*
        ****************************************************************************
        ****************************************************************************
        Create instances of the types of molecular sites
        ****************************************************************************
        ****************************************************************************
        */
         
        SpeciesMethanol species = (SpeciesMethanol)sim.species;
        
        IAtomTypeLeaf type_O  = species.getOxygenType();
        IAtomTypeLeaf type_aC = species.getCarbonType(); 
        IAtomTypeLeaf type_aH = species.getAlphaHydrogenType();
        IAtomTypeLeaf type_H  = species.getHydrogenType();
        IAtomTypeLeaf type_X  = species.getXType();
         
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for calculation of the MOLECULAR potential, U_a_b
        ****************************************************************************
        ****************************************************************************
        */

        
        U_a_b.addPotential(u_O_O,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_O }));
        
        U_a_b.addPotential(u_O_aC,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_aC}));
        U_a_b.addPotential(u_O_aC,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC, type_O }));
        
        U_a_b.addPotential(u_O_aH,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_aH}));
        U_a_b.addPotential(u_O_aH,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH, type_O }));
        
        U_a_b.addPotential(u_O_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_H }));
        U_a_b.addPotential(u_O_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,  type_O }));
        
        
        U_a_b.addPotential(u_aC_aC,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_aC}));
        
        U_a_b.addPotential(u_aC_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_aH}));
        U_a_b.addPotential(u_aC_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_aC}));
        
        U_a_b.addPotential(u_aC_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_H }));
        U_a_b.addPotential(u_aC_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_aC}));
        
        
        U_a_b.addPotential(u_aH_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_aH}));
        
        U_a_b.addPotential(u_aH_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_H }));
        U_a_b.addPotential(u_aH_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_aH}));
        
        U_a_b.addPotential(u_aH_X,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_X }));
        U_a_b.addPotential(u_aH_X,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_X,   type_aH}));
        
        
        U_a_b.addPotential(u_H_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_H }));
        
        
        U_a_b.addPotential(u_X_X,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_X,   type_X }));
        
        
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for simulation, graphics, and data collection
        ****************************************************************************
        ****************************************************************************
        */
        
        sim.integratorOS.setNumSubSteps(1000); // Is this necessary?
        
        
        if (false) { // true to run graphics (and not collect data), false to not run graphics (and collect data)
            sim.box[0].getBoundary().setDimensions(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setDimensions(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            
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
        
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        
        sim.setAccumulatorBlockSize((int)steps);
        
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

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);

        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("chain ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("chain average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("chain overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int nPoints = 2; // number of molecules in simulation (e.g., 2 for B2 calculation)
        public double temperature = 500.0;   // Kelvin
        public long numSteps = 10000;
    }
    
}


