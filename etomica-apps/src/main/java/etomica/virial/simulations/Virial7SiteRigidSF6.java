/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Virial coefficients calculation of SF6
 * 7 sites LJ + rigid without Coulombic interaction 
 * Reference: Samios, Molecular force field investigation for sulfur hexafluoride: A computer simulation study
 * Mayer sampling to evaluate cluster integrals
 * 
 * @author shu
 * 01-18-2013
 */
public class Virial7SiteRigidSF6 {


    public static void main(String[] args) {
    	VirialSF6Param params = new VirialSF6Param();
    	if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
			
		}
        
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;

        System.out.println("SF6 (7 LJ sites without Coulombic interaction, rigid body, from Samios)");
        System.out.println("B"+nPoints+" at T="+temperature+" Kelvin");
		Space space = Space3D.getInstance();
        temperature = Kelvin.UNIT.toSim(temperature);

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
		        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        
        
        // potential for SF6
        double sigmaF = 2.954;
        double sigmaS = 3.246;
        double sigmaSF = (sigmaF +sigmaS ) / 2;// LB rule
        double epsilonF = Kelvin.UNIT.toSim(27.24);
        double epsilonS = Kelvin.UNIT.toSim(163.89);
        double epsilonSF = Math.sqrt(epsilonF * epsilonS);// LB rule
        P2LennardJones p2S = new P2LennardJones(space, sigmaS, epsilonS);
        P2LennardJones p2F = new P2LennardJones(space, sigmaF, epsilonF);
        P2LennardJones pSF = new P2LennardJones(space, sigmaSF , epsilonSF);
        //potential group
        PotentialGroup pGroup = new PotentialGroup(2);
        // fbond and ebond for SF6
        MayerGeneral fTarget= new MayerGeneral(pGroup);
        MayerEGeneral eTarget = new MayerEGeneral(pGroup);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        
        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        
        // species
 //       SpeciesFactory factorySF6 = new SpeciesFactory() {
 //           public ISpecies makeSpecies(ISpace space_) { 
 //           	Species7SiteRigidSF6 species = new Species7SiteRigidSF6(space_);
 //           	species.setConformation(new Conformation7SiteRigidSF6(space_));
 //           	return species;
 //           }
  //      };
    
        // do simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new Species7SiteRigidSF6(space), temperature,refCluster,targetCluster,false);
        sim.box[1].getSampleCluster().value(sim.box[1]);
        sim.integratorOS.setNumSubSteps(1000);
        
        Species7SiteRigidSF6 species = (Species7SiteRigidSF6)sim.getSpecies(0);
        AtomType typeS = species.getSType();
        AtomType typeF = species.getFType();

        pGroup.addPotential(p2S, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeS, typeS}));
        pGroup.addPotential(pSF, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeS, typeF}));
        pGroup.addPotential(p2F, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeF, typeF}));
        pGroup.addPotential(pSF, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeF, typeS}));//switch
        
        
        // graphic part
        if (false) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            
            //set diameters
            DiameterHashByType diameter = new DiameterHashByType(sim); 
            diameter.setDiameter(typeS,1.7);
            diameter.setDiameter(typeF,1.0);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            //set colors
            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.blue);
            colorScheme.setColor(sim.getSpecies(0).getAtomType(1), Color.red);
            
            
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
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
        sim.setAccumulatorBlockSize((int)steps);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        if (true) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(HSB[nPoints]);

 /*       double[] ratioAndError = sim.dvo.getOverlapAverageAndError();
        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
        IData ratioData = ((DataGroup)sim.accumulators[0].getData()).getData(sim.accumulators[0].RATIO.index);
        IData ratioErrorData = ((DataGroup)sim.accumulators[0].getData()).getData(sim.accumulators[0].RATIO_ERROR.index);
        IData averageData = ((DataGroup)sim.accumulators[0].getData()).getData(sim.accumulators[0].AVERAGE.index);
        IData stdevData = ((DataGroup)sim.accumulators[0].getData()).getData(sim.accumulators[0].STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup)sim.accumulators[0].getData()).getData(sim.accumulators[0].ERROR.index);
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[1].RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[1].RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[1].AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[1].STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[1].ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));	
    */
    }


    /**
     * Inner class for parameters
     */
    public static class VirialSF6Param extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 298;
        public long numSteps = 10000;
        public double sigmaHSRef = 5.0;
        public double refFrac = -1;
    }
}                                    
