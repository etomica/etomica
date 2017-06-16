/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import java.awt.Color;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSumNonAdditiveTrimerEnergy;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;

/* 
 * Adapted by Kate from VirialGCPM
 * 
 * Computes only the nonadditive component of either the third, fourth, or fifth virial coefficient for the
 * ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz (JCP 131 064105 2009). 
 * 
 * 
 */


public class VirialCPSHeliumNonAdditive {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        
        double temperatureK; final int nPoints; double sigmaHSRef;
        long steps; int nullRegionMethod;
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperatureK = params.temperature;
            System.out.println(" "+temperatureK);
            steps = params.numSteps;
            sigmaHSRef = params.sigmaHSRef;
            nullRegionMethod = params.nullRegionMethod;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 5) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperatureK = Double.parseDouble(args[1]);
            steps = Integer.parseInt(args[2]);
            sigmaHSRef = Double.parseDouble(args[3]);
            nullRegionMethod = Integer.parseInt(args[4]);
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
        System.out.println("Helium overlap sampling B"+nPoints+"NonAdd at T="+temperatureK+ " K");
        System.out.println("null region method = "+nullRegionMethod);
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        System.out.println(steps+" steps ("+steps/1000+" blocks of 1000)");
        steps /= 1000;

        Space space = Space3D.getInstance();

        
        
        P2HePCKLJS p = new P2HePCKLJS(space);
        P2EffectiveFeynmanHibbs p2 = new P2EffectiveFeynmanHibbs(Space3D.getInstance(), p);
		p2.setTemperature(temperature);		
		p2.setMass(4.002602);
        P3CPSNonAdditiveHe p3NonAdd = new P3CPSNonAdditiveHe(space);
        p3NonAdd.setNullRegionMethod(nullRegionMethod);
    	MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2);
    	ClusterSumNonAdditiveTrimerEnergy targetCluster = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, p3NonAdd, nPoints>3, false);
    	targetCluster.setNo72B2B3NonAdd(true);
    	targetCluster.setTemperature(temperature);

    	
    	MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, null, false);
        refCluster.setTemperature(temperature);


        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), 
                temperature, refCluster,targetCluster, false);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        ///////////////////////////////////////////////
        // Initialize non-overlapped configuration
        ///////////////////////////////////////////////
        
        IAtomList atoms = sim.box[1].getLeafList();
        if (nPoints == 3) {
	        for (int i=1;i<atoms.getAtomCount();i++) {
	        	atoms.getAtom(i).getPosition().setX(0, i*sigmaHSRef);
	        }
        } else if (nPoints == 4) {
	        
	        atoms.getAtom(1).getPosition().setX(0, sigmaHSRef);
	        
	        atoms.getAtom(2).getPosition().setX(0, sigmaHSRef);
	        atoms.getAtom(2).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.getAtom(3).getPosition().setX(1, sigmaHSRef);
	        
        } else if (nPoints == 5) {
        	
        	atoms.getAtom(1).getPosition().setX(0, sigmaHSRef);
        	atoms.getAtom(1).getPosition().setX(1, sigmaHSRef);
        	
	        atoms.getAtom(2).getPosition().setX(0, sigmaHSRef);
	        atoms.getAtom(2).getPosition().setX(1, -sigmaHSRef);
	        
	        atoms.getAtom(3).getPosition().setX(0, -sigmaHSRef);
	        atoms.getAtom(3).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.getAtom(4).getPosition().setX(0, -sigmaHSRef);
	        atoms.getAtom(4).getPosition().setX(1, -sigmaHSRef);
	        
        } else {
        	throw new RuntimeException("Wrong number of points");
        }
        
        
        double r = 4;
        for (int i=1; i<nPoints; i++) {
            Vector v = atoms.getAtom(i).getPosition();
            v.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
//        System.out.println(targetCluster.value(sim.box[1]));
//        System.exit(1);
        
        /*
        IAtomList atoms0 = sim.box[0].getLeafList();
        for (int i=1;i<atoms0.getAtomCount();i++) {
        	atoms0.getAtom(i).getPosition().setX(0, i*sigmaHSRef*1.3);
        } */ 
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
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

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        
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
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 273.15;   // Kelvin
        public long numSteps = 1000000;
        public double sigmaHSRef = 3;
        public int nullRegionMethod = 2;
        public boolean writeRefPref;
    }
}
