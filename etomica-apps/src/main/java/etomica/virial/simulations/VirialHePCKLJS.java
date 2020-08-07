/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * Computes additive virial coefficients using the pair potential for He of Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.
 * 
 * Use QFH boolean to change pair potential to the quadratic Feynman-Hibbs potential.
 * 
 * If determining which option is most efficient via short calculations to estimate standard error, 
 * maintain a 50-50 split of steps between reference and target during data collection with 
 * with adjustStepFreqFrequency = false;
 * 
 * @author kate
 *
 */


public class VirialHePCKLJS {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; int stepsPerBlock; long eqSteps; boolean QFH; boolean adjustStepFrequency;
        if (args.length == 0) {
        	
        	nPoints = params.numMolecules;
            temperatureK = params.temperature;
            steps = params.numSteps;
            stepsPerBlock = params.stepsPerBlock;
            eqSteps = params.eqSteps;
            sigmaHSRef = params.sigmaHSRef;
            QFH = params.QFH;
            adjustStepFrequency = params.adjustStepFrequency;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 8) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperatureK = Double.parseDouble(args[1]);
            steps = Integer.parseInt(args[2]);
            stepsPerBlock = Integer.parseInt(args[3]);
            eqSteps = Integer.parseInt(args[4]);
            sigmaHSRef = Double.parseDouble(args[5]);
            QFH = Boolean.parseBoolean(args[6]);
            adjustStepFrequency = Boolean.parseBoolean(args[7]);
            params.writeRefPref = true;
        	
        }  else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
        

        int numSubSteps = 1000;

        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        

        System.out.println("Target diagram: B"+nPoints+" for helium pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        if (QFH) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        //Next line not needed because energy in Kelvin
        //temperature = Kelvin.UNIT.toSim(temperature);

        
        Space space = Space3D.getInstance();

        MayerGeneralSpherical fTarget; MayerESpherical eTarget;
        if (QFH) {
        	P2HePCKLJS p20 = new P2HePCKLJS(space);
        	P2EffectiveFeynmanHibbs p2 = new P2EffectiveFeynmanHibbs(Space3D.getInstance(),p20);
			p2.setTemperature(temperature);		
			p2.setMass(4.002602);
			fTarget = new MayerGeneralSpherical(p2);
	    	eTarget = new MayerESpherical(p2);
        } else {
        	P2HePCKLJS p2 = new P2HePCKLJS(space);
        	fTarget = new MayerGeneralSpherical(p2);
        	eTarget = new MayerESpherical(p2);
        }
    	
    	ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, false);
   	    ClusterWeight sampleCluster1 = ClusterWeightAbs.makeWeightCluster(targetCluster);

    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        ClusterWeight refSample = ClusterWeightAbs.makeWeightCluster(refCluster);

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        sampleCluster1.setTemperature(temperature);
        refSample.setTemperature(temperature);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), 
                temperature, refCluster,targetCluster, false);
        
        System.out.println((steps*1000)+" steps ("+steps+" IntegratorOverlap steps of 1000 steps)");
        sim.integratorOS.setNumSubSteps(1000);
        
        // Keep 50-50 split between reference and target
        sim.integratorOS.setAdjustStepFraction(adjustStepFrequency);
       
        sim.setAccumulatorBlockSize(stepsPerBlock);
        System.out.println(stepsPerBlock+" steps per block");
        
        System.out.println();

        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesSpheresMono species = (SpeciesSpheresMono)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0, false);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, eqSteps); // 5000 IntegratorOverlap steps = 5e6 steps
        System.out.println((eqSteps*1000) + " equilibration steps (" + eqSteps + " Integrator Overlap Steps)"); 
        
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);

        long t1 = System.currentTimeMillis();

        sim.getController().actionPerformed();
        System.out.println();
        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        System.out.println();
        double[] ratioAndError = sim.dvo.getAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
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
        
        System.out.println();
        System.out.println("cm"+((nPoints-1)*3)+"/mol"+(nPoints-1)+": ");
        System.out.println("abs average: "+ratio*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1)+", error: "+error*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {

        public int numMolecules = 4;
        public double temperature = 250;
        public long numSteps = 1000;
        public int stepsPerBlock = 1000;
        public long eqSteps=1000;
        public double sigmaHSRef = 3;
        public boolean writeRefPref = false;
        public boolean QFH = true;
        public boolean adjustStepFrequency = true;
    }
}
