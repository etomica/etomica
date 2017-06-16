/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * Anthracene model , pure 
 * support for CO2/An mixture calculation
 * 3 site model from Iwai, this is 464 model
 * Mayer sampling to evaluate cluster integrals
 * @author shu
 * 
 * March,9,2011
 */
public class VirialAnthracene464 {


    public static void main(String[] args) {
        VirialAnthracene545Param params = new VirialAnthracene545Param();
        if (args.length > 0) {
            ReadParameters readParameters = new ReadParameters(args[0], params);
            readParameters.readParameters();
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;

        System.out.println("Anthracene(464) overlap sampling B"+nPoints+" at T="+temperature+" Kelvin");
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
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        
        
        // potential for Anthracene
        double sigmaC = 5.502;
        double sigmaCH = 5.502;
        double sigmaCCH = 5.502;
        double epsilonC = Kelvin.UNIT.toSim(360.266356);
        double epsilonCH = Kelvin.UNIT.toSim(145.5712827);
        double epsilonCCH = Math.sqrt(epsilonC*epsilonCH);
        P2LennardJones p2C = new P2LennardJones(space, sigmaC, epsilonC);
        P2LennardJones p2CH = new P2LennardJones(space, sigmaCH, epsilonCH);
        P2LennardJones pCCH = new P2LennardJones(space,sigmaCCH , epsilonCCH);
        //potential group
        PotentialGroup pTarget = new PotentialGroup(2);
        // f-bond and e-bond for anthracene
        MayerGeneral fTarget= new MayerGeneral(pTarget);
        MayerEGeneral eTarget = new MayerEGeneral(pTarget);
        
        // cluster
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
    
        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        
        // species anthracene
        SpeciesFactory factoryAn = new SpeciesFactory() {
            public ISpecies makeSpecies(Space space) {
            	SpeciesAnthracene3site464 species = new SpeciesAnthracene3site464(space);
                      return species;
            }
        };
    
    // do simulation
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,factoryAn, temperature,refCluster,targetCluster,false);
        sim.box[1].getSampleCluster().value(sim.box[1]);
        sim.integratorOS.setNumSubSteps(1000);
                
        //put the species in the box
        SpeciesAnthracene3site464 speciesAn = (SpeciesAnthracene3site464)sim.getSpecies(0);
        sim.integratorOS.setNumSubSteps(1000);

        AtomType typeC = speciesAn.getCType();
        AtomType typeCH = speciesAn.getCHType();

        pTarget.addPotential(p2C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));
        pTarget.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH}));
        pTarget.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeCH}));
        pTarget.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeC}));//switch
        // graphic part
        if (false) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
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
            sim.integratorOS.setStepFreq0(refFrac);
            sim.integratorOS.setAdjustStepFreq(false);
        }

        if (true) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
        IData ratioData = ((DataGroup) sim.accumulators[0].getData()).getData(AccumulatorRatioAverageCovariance.RATIO.index);
        IData ratioErrorData = ((DataGroup) sim.accumulators[0].getData()).getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index);
        IData averageData = ((DataGroup) sim.accumulators[0].getData()).getData(AccumulatorAverage.AVERAGE.index);
        IData stdevData = ((DataGroup) sim.accumulators[0].getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup) sim.accumulators[0].getData()).getData(AccumulatorAverage.ERROR.index);
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));

        ratioData = ((DataGroup) sim.accumulators[1].getData()).getData(AccumulatorRatioAverageCovariance.RATIO.index);
        ratioErrorData = ((DataGroup) sim.accumulators[1].getData()).getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index);
        averageData = ((DataGroup) sim.accumulators[1].getData()).getData(AccumulatorAverage.AVERAGE.index);
        stdevData = ((DataGroup) sim.accumulators[1].getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        errorData = ((DataGroup) sim.accumulators[1].getData()).getData(AccumulatorAverage.ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
	}

    /**6
     * Inner class for parameters
     */
    public static class VirialAnthracene545Param extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 400;
        public long numSteps = 100000;
        public double sigmaHSRef = 5*1.40;
        public double refFrac = -1;
    }
}                                    
