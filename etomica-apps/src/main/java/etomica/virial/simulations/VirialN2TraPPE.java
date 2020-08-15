/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Nitrogen;
import etomica.config.IConformation;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.P2CO2EMP;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * TraPPE N2 model 
 * Potoff and Siepmann. Vapor-Liquid Eq of Mixtures Containing Alkanes, Carbon Dioxide and Nitrogen, 1999
 * N site:LJ+ negative partial charge, COM:positive charge
 * Mayer sampling to evaluate cluster integrals
 * 
 * @author shu
 * 01-18-2013
 */
public class VirialN2TraPPE {

	public static void main(String[] args) {
       
		VirialN2Param params = new VirialN2Param();
		if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
			
		}
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        
        double sigmaA = params.sigmaA;
        double sigmaN = params.sigmaN;
        double epsilonN = params.epsilonN;
        double sigmaAN = (sigmaA + sigmaN ) / 2.0;//LB rule
        double charge = params.charge;
        System.out.println("N2(TraPPE) overlap sampling B"+nPoints+" at T="+temperature+" Kelvin");
        temperature = Kelvin.UNIT.toSim(temperature);

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
//        double beta = 107.8 * 0.5 * Math.PI/ 180.0;
//        double asin =  Math.asin(    Math.sin(beta ) * 2 / Math.sqrt(3)    );
//        double a = Math.PI - asin; 
//        double real  = a * 180.0 / Math.PI;
//        double real1 = (Math.PI - Math.asin(    Math.sin(107.8 * 0.5 * Math.PI/ 180.0 ) * 2 / Math.sqrt(3)    )) * 180.0 / Math.PI;
  //      System.out.println("real is "+real);
 //       System.out.println("beta is "+beta);
  //      System.out.println("real1 is "+real1);
        
   //     double beta1 = Math.PI- 110.7 * Math.PI / 180.0;
   //     double real2 = Math.acos(   1 - 1.5 * Math.sin(beta1)* Math.sin(beta1));
   //     System.out.println("real2 is "+real2 * 180.0 / Math.PI);

        
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println(" twice B2HS:"+  2 * HSB[2]);
        System.out.println(" three times B3HS:"+  3 * HSB[3]);
        System.out.println("B4HS:"+  HSB[4]);
        System.out.println("four times the B4HS:"+ 4* HSB[4]);
        System.out.println("five times the B5HS:"+ 5* HSB[5]);
        System.out.println("B6HS:"+  HSB[6]);
        
        Space space = Space3D.getInstance();
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);        
        // P2CO2EMP(ISpace space, double sigmaC, double sigmaCO, double sigmaO, double epsilonC, double epsilonCO, double epsilonO, double chargeC) 
        P2CO2EMP pN2 = new P2CO2EMP(space, sigmaA, sigmaAN, sigmaN, 0.0, 0.0, Kelvin.UNIT.toSim(epsilonN), Electron.UNIT.toSim(charge));
        MayerGeneral fN2 = new MayerGeneral(pN2);
        MayerEGeneral eN2 = new MayerEGeneral(pN2);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fN2, nPoints>3, eN2, true);
        
        refCluster.setTemperature(temperature);
        targetCluster.setTemperature(temperature);
        
        // nitrogen conformation
        final IConformation conformation = new IConformation() {
        	public void initializePositions(IAtomList atomList) {
        		// atoms are N-COM("A")-N, 1-0-2 
        		double bondL = 1.10;// length of N--N
                atomList.get(0).getPosition().E(0);// COM("A")
                atomList.get(1).getPosition().E(0);// Nitrogen atom on the left
                atomList.get(1).getPosition().setX(0, -0.5 * bondL);
                atomList.get(2).getPosition().E(0);// Nitrogen atom on the right
                atomList.get(2).getPosition().setX(0, +0.5 * bondL);
            }
        };
        
        SpeciesSpheresHetero species = new SpeciesSpheresHetero(space, new IElement[]{new ElementSimple("A"), Nitrogen.INSTANCE});
        species.setChildCount(new int[]{1,2});
        species.setConformation(conformation);
        
        //simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species, temperature,refCluster,targetCluster, false);
//        sim.box[1].getSampleCluster().value(sim.box[1]);
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;
	// graphical part
        if (false) {
            double size = 5.0;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            //set diameters
            DiameterHashByType diameter = new DiameterHashByType();
            diameter.setDiameter(species.getAtomType(0), 0.2);
            diameter.setDiameter(species.getAtomType(1), 0.2);


            simGraphic.getDisplayBox(sim.box[0]).setDiameterHash(diameter);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);

            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.gray);
            colorScheme.setColor(sim.getSpecies(0).getAtomType(1), Color.cyan);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController2().addActivity(new ActivityIntegrate2(sim.integratorOS));
            if (Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0) {
                throw new RuntimeException("Oops");
            }

            return;
        }
        
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }

        System.out.println("equilibration finished");

        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.integratorOS.getMoveManager().setEquilibrating(false);

        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }
        if (refFrac >= 0) {
        	 sim.integratorOS.setRefStepFraction(refFrac);
             sim.integratorOS.setAdjustStepFraction(false);
        }

        if (false) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.getController2().getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), 1000);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(HSB[nPoints]);
        /*
        // convert results' units
        Unit convertFromMacro = new UnitRatio(Liter.UNIT, Mole.UNIT);
        double ref = 21.3; 
        System.out.println("reference B2 at T =600K :" + convertFromMacro.toSim(ref/1000));
        Unit convertFromSimUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);       
        System.out.println("convert simulation units to macro units:");

        double simResult_B2= -94.54719724658284  ;
        System.out.println("B2 ********************** (A3/atom)  :  " + simResult_B2);
        System.out.println("convert  to L/mol :  " + convertFromSimUnit.toSim(simResult_B2));
        System.out.println( "or in  L/mol :  " + simResult_B2 * 6.02214129/ 10000.0);
        System.out.println( "or in  cm3/mol :  " + simResult_B2 * 6.02214129 / 10.0);
        
        double simResult_B3 = 3608.75979282736     ;
        System.out.println("my B3 (A3/atom)^2  :  " + simResult_B3);
        System.out.println( "convert  to L2/mol2 :  "+ simResult_B3 * 6.02214129 * 6.02214129 / 10000.0 / 10000.0);/////?????
        System.out.println( "convert  to cm6/mol2 :  "+ simResult_B3 * 6.02214129 * 6.02214129 / 100.0);/////?????
        
        Unit fromMacro_epsilon = new UnitRatio(Joule.UNIT, Mole.UNIT);
        System.out.println("in sim unit the epsilon is : " + fromMacro_epsilon.toSim(0.60807 * 1000));
        System.out.println("in sim unit the epsilon is : " + fromMacro_epsilon.toSim(0.586 * 1000));

        Unit fromSim_epsilon = new UnitRatio(Mole.UNIT , Joule.UNIT);
        System.out.println("***epsilon is : " + fromSim_epsilon.toSim(70.60));
 //       System.out.println(" epsilon is : " + fromMacro_epsilon.toSim(0.60807 * 1000));
   //     System.out.println(" is : " + fromMacro_epsilon.toSim(0.60807 * 1000));
  //      System.out.println(" is : " + fromMacro_epsilon.toSim(0.60807 * 1000));
 //       System.out.println("in sim unit the epsilon is : " + fromMacro_epsilon.toSim(0.60807 * 1000));
*/
	}

    /**
     * Inner class for parameters
     */
    public static class VirialN2Param extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 291.41;
        public long numSteps = 1000000000;
        public double sigmaHSRef = 5.0;
        public double sigmaA = 1.0;// set arbitrarily
        public double sigmaN = 3.31;
        public double epsilonN = 36.0;
        public double charge = 0.964;// in the center of mass
        public double refFrac = -1;
    }
    
 
    
}
