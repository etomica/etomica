/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.AlkaneEH.SpeciesMethane;
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
 * TraPPE-EH methane, CH4  
 * Mayer sampling to evaluate cluster integrals
 * 
 * @author shu
 * 02-12-2013
 */
public class VirialMethaneEH {


    public static void main(String[] args) {
        VirialMethaneParam params = new VirialMethaneParam();
        if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
		}
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;

        System.out.println("new , CH4 (TraPPE-EH) overlap sampling B"+nPoints+" at T="+temperature+" Kelvin");
        System.out.println("new angle! ");
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
     
        Space space = Space3D.getInstance();
        // ref system
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        double sigmaH = 3.31;// "middle point of CH bond"
        double sigmaC = 3.31;
        double sigmaCH = (sigmaH+sigmaC)/2;
        double epsilonH = Kelvin.UNIT.toSim(15.3);
        double epsilonC = Kelvin.UNIT.toSim(0.01);
        double epsilonCH = Math.sqrt((epsilonH * epsilonC ));
        //System.out.println(epsilonH*16+epsilonC+epsilonCH*8);
        P2LennardJones p2H = new P2LennardJones(space, sigmaH, epsilonH);// H-H
        P2LennardJones p2C = new P2LennardJones(space, sigmaC, epsilonC);//C-C
        P2LennardJones p2CH = new P2LennardJones(space, sigmaCH, epsilonCH);//C-H
        PotentialGroup pTargetGroup = new PotentialGroup(2);
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup);
        
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        SpeciesMethane species = new SpeciesMethane(space);

        //simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species, temperature,refCluster,targetCluster,false);

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");

        steps /= 1000;

        AtomType typeC = species.getAtomType(0);
        AtomType typeH = species.getAtomType(1);

        // build methane potential
        pTargetGroup.addPotential(p2C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));//C-C
        pTargetGroup.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeH}));//C-H
        pTargetGroup.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeC}));//H-C
        pTargetGroup.addPotential(p2H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));//H-H

        if (false) {
            double size = 5.0;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            //set diameters
            DiameterHashByType diameter = new DiameterHashByType(sim); 
            diameter.setDiameter(species.getAtomType(0),0.4);
            diameter.setDiameter(species.getAtomType(1),0.2);
            
            simGraphic.getDisplayBox(sim.box[0]).setDiameterHash(diameter);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            
            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.gray);
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
        
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
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

      
	}

    /**
     * Inner class for parameters
     */
    public static class VirialMethaneParam extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 300;
        public long numSteps = 10000000;
        public double sigmaHSRef = 6;
        public double refFrac = -1;
    }
}
