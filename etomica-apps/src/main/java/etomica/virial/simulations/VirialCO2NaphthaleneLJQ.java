/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2LJQQ;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.CompoundUnit;
import etomica.units.Coulomb;
import etomica.units.Kelvin;
import etomica.units.Meter;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * single site for CO2 and single site for naphthalene
 * mixture of CO2 and Naphthalene
 * 
 * the code is modified from virialco2SKSmix and virialCO2NaphthaleneTraPPE
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 * parameters are from albo et al
 * @author shu
 * Dec.15,2010
 * 
 */
public class VirialCO2NaphthaleneLJQ {

    public static void main(String[] args) {
    	VirialCO2NaphthaleneLJQParam params = new VirialCO2NaphthaleneLJQParam();
        if (args.length > 0) {
            ReadParameters readParameters = new ReadParameters(args[0], params);
            readParameters.readParameters();
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int[] nTypes = params.nTypes;
        double refFrac = params.refFrac;
        int sum = 0;
        for (int i=0; i<nTypes.length; i++) {
            if (nTypes[i] == 0) {// pure 
                throw new RuntimeException("for pure component, use a different class");
            }
            sum += nTypes[i];
        }
        if (sum != nPoints) {
            throw new RuntimeException("Number of each type needs to add up to nPoints");
        }
        
        System.out.println("CO2 + Naphthalene LJQ overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
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
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        
      //Target system CO2  potential from P2LJQQ
        double sigmaCO2 = 3.748 ;
        // double epsilon = 215 ; 
        double epsilonCO2 = Kelvin.UNIT.toSim(215);
     // q :Coulombs * m^2, Q : sim unit
        double QCO2 = new CompoundUnit(new Unit[]{Coulomb.UNIT,Meter.UNIT}, new double [] {1,2}).toSim(-1.367e-39);
        double QCO22 = QCO2 * QCO2;
        P2LJQQ pCO2 = new P2LJQQ(space, sigmaCO2, epsilonCO2, QCO22);
        //CO2 potential need to know the temperature!
        pCO2.setTemperature(temperature);
                //mayer function , f-bond and e-bond
        MayerGeneralSpherical fCO2 = new MayerGeneralSpherical(pCO2);
        MayerESpherical eCO2 = new MayerESpherical(pCO2);
           
        
        //this is the Na-Na potential from P2LJQQ
        double sigmaNa = 5.5806 ;
        // double epsilon = 530 ; 
        double epsilonNa = Kelvin.UNIT.toSim(530);
     // q :Coulombs * m^2, Q : sim unit
        double QNa = new CompoundUnit(new Unit[]{Coulomb.UNIT,Meter.UNIT}, new double [] {1,2}).toSim(4.503e-39);
        double QNa2 = QNa * QNa;
        P2LJQQ pNa = new P2LJQQ(space, sigmaNa, epsilonNa, QNa2);
        //potential needs to know T
        pNa.setTemperature(temperature);
        //mayer function , f-bond and e-bond
        MayerGeneralSpherical fNa = new MayerGeneralSpherical(pNa);
        MayerESpherical eNa = new MayerESpherical(pNa);

        // this is CO2-Na interaction potential, the 1st is CO2, 2nd is Na

        double sigmaCO2Na =  (sigmaCO2 + sigmaNa)/2;
        double epsilonCO2Na = Math.sqrt(epsilonCO2*epsilonNa);
        double QCO2Na = QCO2 * QNa; 
        P2LJQQ pCO2Na = new P2LJQQ(space, sigmaCO2Na, epsilonCO2Na, QCO2Na);
        pCO2Na.setTemperature(temperature);

        MayerGeneralSpherical fCO2Na= new MayerGeneralSpherical(pCO2Na);
        MayerESpherical eCO2Na = new MayerESpherical(pCO2Na);
        
        
        // define target cluster, here we use mixture's cluster   i am not sure about this !
        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fCO2,fCO2Na},{fCO2Na,fNa}},
        new MayerFunction[][]{{eCO2,eCO2Na},{eCO2Na,eNa}}, nTypes);
        targetCluster.setTemperature(temperature);
        

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");

        


        // now is the simulation!!!
        //SpeciesFactorySpheres factoryCO2 = new SpeciesFactorySpheres();
        //SpeciesFactorySpheres factoryNa = new SpeciesFactorySpheres();
      final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactorySpheres(),temperature,new ClusterAbstract[]{refCluster,targetCluster},
              new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},false);
        
   
        //graphic part
        if (false) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
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
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
        if (refFrac >= 0) {
            sim.integratorOS.setStepFreq0(refFrac);
            sim.integratorOS.setAdjustStepFreq(false);
        }
        
        

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
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
        IData ratioData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverageCovariance.RATIO.index);
        IData ratioErrorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index);
        IData averageData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorAverage.AVERAGE.index);
        IData stdevData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorAverage.ERROR.index);
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverageCovariance.RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorAverage.AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorAverage.ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
	}

    /**
     * Inner class for parameters
     */
    public static class VirialCO2NaphthaleneLJQParam extends ParameterBase {
        public int nPoints = 3;
        public double temperature = 328.15;
        public long numSteps = 10000;
        public double sigmaHSRef = 7;
       //composition
        public int[] nTypes = new int[]{2,1};
        public double refFrac = -1;
    }
}
