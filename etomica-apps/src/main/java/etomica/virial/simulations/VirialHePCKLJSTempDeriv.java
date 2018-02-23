/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import java.awt.Color;

import etomica.atom.IAtomList;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.MayerD2FDT2Spherical;
import etomica.virial.MayerDFDTSpherical;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Computes first or second temperatures (see firstDerivative variable) of additive virial coefficients using the 
 * ab initio pair potential for helium from Przybytek et al. (2010) (Phys. Rev. Lett. 104, 183003). 
 * 
 * Use the boolean QFH to select whether the quadratic Feynman-Hibbs modification to the potential is used.
 * 
 * @author kate 
 *
 */


public class VirialHePCKLJSTempDeriv {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; boolean firstDerivative; boolean QFH;
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperatureK = params.temperature;
            steps = params.numSteps;
            sigmaHSRef = params.sigmaHSRef;
            firstDerivative = params.firstDerivative;
            QFH = params.QFH;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 6) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperatureK = Double.parseDouble(args[1]);
            steps = Integer.parseInt(args[2]);
            sigmaHSRef = Double.parseDouble(args[3]);
            firstDerivative = Boolean.parseBoolean(args[4]);
            QFH = Boolean.parseBoolean(args[5]);
            params.writeRefPref = true;
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
        

        int numSubSteps = 1000;

        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        if (firstDerivative) {
        	System.out.println("Helium overlap sampling dB"+nPoints+"Add/dT at T="+temperatureK+ " K");
        } else {
        	System.out.println("Helium overlap sampling d2B"+nPoints+"Add/dT2 at T="+temperatureK+ " K");
        }
        if (QFH) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        //Next line not needed because energy in Kelvin
        //temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println(steps+" steps ("+steps/1000+" blocks of 1000)");
        steps /= 1000;

        Space space = Space3D.getInstance();

        MayerGeneralSpherical fTarget; MayerESpherical eTarget; 
        MayerDFDTSpherical dfdTTarget; MayerD2FDT2Spherical df2dT2Target;
        if (QFH) {
        	P2HePCKLJS p20 = new P2HePCKLJS(space);
        	P2EffectiveFeynmanHibbs p2 = new P2EffectiveFeynmanHibbs(Space3D.getInstance(),p20);
			p2.setTemperature(temperature);		
			p2.setMass(4.002602);
			fTarget = new MayerGeneralSpherical(p2);
	    	eTarget = new MayerESpherical(p2);
	    	dfdTTarget = new MayerDFDTSpherical(p2);
        	df2dT2Target = new MayerD2FDT2Spherical(p2);
        } else {
        	P2HePCKLJS p2 = new P2HePCKLJS(space);
        	fTarget = new MayerGeneralSpherical(p2);
        	eTarget = new MayerESpherical(p2);
        	dfdTTarget = new MayerDFDTSpherical(p2);
        	df2dT2Target = new MayerD2FDT2Spherical(p2);
        }
        


    	
    	VirialDiagrams targetCluster1 = new VirialDiagrams(nPoints, false, false);
    	ClusterSum targetCluster;
    	if (firstDerivative) {
    		targetCluster = targetCluster1.makeVirialClusterTempDeriv(fTarget, eTarget, dfdTTarget);
    	} else {
    		targetCluster = targetCluster1.makeVirialClusterSecondTemperatureDerivative(fTarget, eTarget, dfdTTarget, df2dT2Target);
    	}
    	
        
        
    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
       
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
       
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactorySpheres(), 
                temperature,refCluster,targetCluster, false);
        
        IAtomList atoms = sim.box[1].getLeafList();
        if (QFH) {
        	
        } else {
        	for (int i = 1; i<atoms.size(); i++) {
            	atoms.get(i).getPosition().setX(0, i*10);
            }
        }
        
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
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
//                    sim.initRefPref(null, 0);
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
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());
        
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

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
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
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[0].RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[0].RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[0].AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[0].STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(sim.accumulators[0].ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        System.out.println();
        System.out.println("cm"+((nPoints-1)*3)+"/mol"+(nPoints-1)+"/K: ");
        
        double convertBeta = temperature/temperatureK;
        if (!firstDerivative) {
        	convertBeta*=convertBeta;
        }
        System.out.println("abs average: "+ratio*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1)*convertBeta+", error: "+error*HSB[nPoints]*Math.pow(Constants.AVOGADRO*1e-24,nPoints-1)*convertBeta);
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 50.0;   // Kelvin
        public long numSteps = 100000;
        public double sigmaHSRef = 3;
        public boolean firstDerivative = true;
        public boolean QFH = false;
        public boolean writeRefPref;
    }
}
