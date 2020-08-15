/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.potential.P2ArgonAziz1993;
import etomica.potential.P2QChem;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.virial.B2ForSphericallySymmetricUByTrapezoidRule;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationClusterMove;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * This class used the direct-sampling formulation of Mayer-Sampling Monte Carlo to compute the second
 * virial coefficient of argon, as described with a quantum potential (via Q-Chem), at a particular temperature.
 * The argon.in file (written via P2QChem.java) controls which quantum-mechanical description is employed.
 *
 * The pair (true pair) potential of Aziz (1993) is employed as the reference system.  
 *  
 * To test this simulation within Eclipse, one can set the boolean variable "practice" to true.  This causes the
 * pair potential of Aziz to be used for both the reference and target systems.  Q-Chem can only be used on CCR's U2.
 * 
 * Kate Shaul, 2009/2010
 */
public class DirectSamplingTargetQCArReferenceAzizAr {


    public static void main(String[] args) {
    	
    	// ***********************************************************************************************************
        // If (practice), then the Aziz (1993) potential is used for both the target and reference potential.
        // ***********************************************************************************************************
    	boolean practice = false;

    	VirialLJParam params = new VirialLJParam();

        double temperature; final int nPoints; int log10Steps;
        
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperature = params.temperature;
            log10Steps = params.log10Steps;
            
        } else if (args.length == 3) {

        	nPoints = Integer.parseInt(args[0]);
        	temperature = Double.parseDouble(args[1]);
            log10Steps = Integer.parseInt(args[2]);      
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments.");
        }

 
        long steps = (long)Math.pow(10, log10Steps);

        
        // ***********************************************************************************************************
        // Computation of the reference B2 value
        // ***********************************************************************************************************
        Space space = Space3D.getInstance();
        P2ArgonAziz1993 p2 = new P2ArgonAziz1993(space);
        B2ForSphericallySymmetricUByTrapezoidRule b2ref = new B2ForSphericallySymmetricUByTrapezoidRule();
        double[] results = b2ref.getResults(p2, 10, 100, temperature, false);
        double ref = results[0];
        double refError = results[1];
        
        System.out.println("Argon B" + nPoints);
        System.out.println("T = " + temperature );
    
        System.out.println();
        System.out.println("Direct sampling with pair potential of Aziz (1993) as reference");
        System.out.println("Reference B" + nPoints + " (T = " +temperature +" K) = "+ ref + " +/- " + refError);
        System.out.println();
        System.out.println("10^" + log10Steps + " steps");
        System.out.println();

        P2ArgonAziz1993 p2Ref = new P2ArgonAziz1993(space); // computes "energy" in Kelvin 
        
        MayerGeneralSpherical fRef = new MayerGeneralSpherical(p2Ref);
        MayerESpherical eRef = new MayerESpherical(p2Ref);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        
        ClusterAbstract[] targetCluster = new ClusterAbstract[1];
        
        
        if (practice) {
        	// Q-Chem (and thus P2QChem) cannot be used outside of CCR
        	MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2Ref);
            MayerESpherical eTarget = new  MayerESpherical(p2Ref);
            targetCluster[0]= Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        } else {
        	P2QChem p2QChem = new P2QChem(space); // computes "energy" in Kelvin 
        	MayerGeneral fTarget = new MayerGeneral(p2QChem);
        	MayerEGeneral eTarget = new  MayerEGeneral(p2QChem);
            targetCluster[0]= Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        }
        
        targetCluster[0].setTemperature(temperature);

        // The sampleCluster is the system used to sample configuration space. 
        ClusterWeight sampleCluster = ClusterWeightAbs.makeWeightCluster(refCluster);
        
        SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
		
        final SimulationVirial sim = new SimulationVirial(space, speciesFactory, temperature, sampleCluster,refCluster,targetCluster);
        
        ConfigurationClusterMove clusterMove = new ConfigurationClusterMove(space, sim.getRandom());
        clusterMove.initializeCoordinates(sim.box);
        
        sim.equilibrate(steps/40);
        
        System.out.println("Equilibration finished.");
        
        System.out.println();
        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize());
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), steps);

//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        
        System.out.println();
        System.out.println("Reference average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
                           +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
                           +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);
      
        System.out.println();
    		
        double	ratio = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];
        double  error = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];
        
    	System.out.println("Target average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1]
                               +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[1]
                               +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1]);
        
    	System.out.println();
            
        System.out.println("Ratio average: "+ratio+", error: "+error);
        
        System.out.println();

        // ***********************************************************************************************************
        // Error from reference B2 value is neglected.
        // ***********************************************************************************************************
        System.out.println("B"+nPoints+": "+ratio*ref+", error: "+ error*ref);
        
        System.out.println();

	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 200; // Kelvin
        public int log10Steps = 5; // steps = (long)Math.pow(10, 6);
        
    }
}
