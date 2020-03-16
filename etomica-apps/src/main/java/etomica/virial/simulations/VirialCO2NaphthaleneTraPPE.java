/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.models.co2.SpeciesTraPPECO2;
import etomica.potential.P2CO2TraPPE;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * TraPPE
 * mixture of CO2 and Naphthalene
 * 
 * the code is modified from virialco2SKSmix
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 * 
 * @author shu
 * 10.22.2010
 * 
 */
public class VirialCO2NaphthaleneTraPPE {

    public static void main(String[] args) {
        VirialCO2TraPPEParam params = new VirialCO2TraPPEParam();
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
        
        System.out.println("CO2 + Naphthalene overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
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
        
    // CO2 TraPPE potential
        P2CO2TraPPE pCO2 = new P2CO2TraPPE(space);
        MayerGeneral fCO2 = new MayerGeneral(pCO2);
        MayerEGeneral eCO2 = new MayerEGeneral(pCO2);
    
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        
        //this is the Na-Na potential 
        double sigmaC = 3.70;
        double sigmaCH = 3.695;
        double sigmaCCH = (sigmaC +sigmaCH ) / 2;
        double epsilonC = Kelvin.UNIT.toSim(30.0);
        double epsilonCH = Kelvin.UNIT.toSim(50.5);
        double epsilonCCH = Math.sqrt(epsilonC*epsilonCH);
        P2LennardJones p2C = new P2LennardJones(space, sigmaC, epsilonC);
        P2LennardJones p2CH = new P2LennardJones(space, sigmaCH, epsilonCH);
        P2LennardJones pCCH = new P2LennardJones(space,sigmaCCH , epsilonCCH);
       
        PotentialGroup pNa = new PotentialGroup(2);
        MayerGeneral fNa= new MayerGeneral(pNa);
        MayerEGeneral eNa = new MayerEGeneral(pNa);
        
          // this is CO2-Na interaction potential, the 1st is CO2, 2nd is Na
        P2LennardJones pC_C = new P2LennardJones(space, 0.5*(sigmaC+pCO2.getSigmaC()), Math.sqrt(epsilonC*pCO2.getEpsilonC()));
        P2LennardJones pC_CH = new P2LennardJones(space, 0.5*(sigmaCH+pCO2.getSigmaC()), Math.sqrt(epsilonCH*pCO2.getEpsilonC()));
        P2LennardJones pO_C = new P2LennardJones(space, 0.5*(sigmaC+pCO2.getSigmaO()), Math.sqrt(epsilonC*pCO2.getEpsilonO()));
        P2LennardJones pO_CH = new P2LennardJones(space, 0.5*(sigmaCH+pCO2.getSigmaO()), Math.sqrt(epsilonCH*pCO2.getEpsilonO()));
   
        PotentialGroup pCO2Na = new PotentialGroup(2);
        MayerGeneral fCO2Na= new MayerGeneral(pCO2Na);
        MayerEGeneral eCO2Na = new MayerEGeneral(pCO2Na);
        
        
        // define target cluster, here we use mixture's cluster
        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fCO2,fCO2Na},{fCO2Na,fNa}},
                new MayerFunction[][]{{eCO2,eCO2Na},{eCO2Na,eNa}}, nTypes);
        targetCluster.setTemperature(temperature);
        

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		// initialize
      
      
        // put SpeciesFactoryCO2 and SpeciesFactoryNa here and these classes contain conformations already
        
        // this is CO2
        SpeciesFactory factoryCO2 = new SpeciesFactory() {
            public ISpecies makeSpecies(Space space) { //declare
            	SpeciesTraPPECO2 species = new SpeciesTraPPECO2(space);
                      return species;
            }
        };

        // this is for Na
        SpeciesFactory factoryNa = new SpeciesFactory() {
            public ISpecies makeSpecies(Space space) {
            	SpeciesTraPPENaphthalene species = new SpeciesTraPPENaphthalene(space);
                      return species;
            }
        };
        
        
          // now is the simulation!!!
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactory[]{factoryCO2,factoryNa}, nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},true);
        
        //put the species in the box
        SpeciesTraPPENaphthalene speciesNa = (SpeciesTraPPENaphthalene)sim.getSpecies(1);
        SpeciesTraPPECO2 speciesCO2 = (SpeciesTraPPECO2)sim.getSpecies(0);
                sim.integratorOS.setNumSubSteps(1000);

        AtomType typeC = speciesNa.getCType();
        AtomType typeCH = speciesNa.getCHType();
        AtomType typeC_CO2 = speciesCO2.getCarbonType();
        AtomType typeO = speciesCO2.getOxygenType();

        pNa.addPotential(p2C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));
        pNa.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH}));
        pNa.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeCH}));
        pNa.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeC}));//switch
        pCO2Na.addPotential(pC_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CO2, typeC}));
        pCO2Na.addPotential(pC_CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CO2, typeCH}));
        pCO2Na.addPotential(pO_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeC}));
        pCO2Na.addPotential(pO_CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeCH}));

        
        if (false) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
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
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        if (refFrac >= 0) {
            sim.integratorOS.setStepFreq0(refFrac);
            sim.integratorOS.setAdjustStepFreq(false);
        }

        if (false) {
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
        IData ratioData = ((DataGroup) sim.accumulators[0].getData()).getData(sim.accumulators[0].RATIO.index);
        IData ratioErrorData = ((DataGroup) sim.accumulators[0].getData()).getData(sim.accumulators[0].RATIO_ERROR.index);
        IData averageData = ((DataGroup) sim.accumulators[0].getData()).getData(sim.accumulators[0].AVERAGE.index);
        IData stdevData = ((DataGroup) sim.accumulators[0].getData()).getData(sim.accumulators[0].STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup) sim.accumulators[0].getData()).getData(sim.accumulators[0].ERROR.index);
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));

        ratioData = ((DataGroup) sim.accumulators[1].getData()).getData(sim.accumulators[1].RATIO.index);
        ratioErrorData = ((DataGroup) sim.accumulators[1].getData()).getData(sim.accumulators[1].RATIO_ERROR.index);
        averageData = ((DataGroup) sim.accumulators[1].getData()).getData(sim.accumulators[1].AVERAGE.index);
        stdevData = ((DataGroup) sim.accumulators[1].getData()).getData(sim.accumulators[1].STANDARD_DEVIATION.index);
        errorData = ((DataGroup) sim.accumulators[1].getData()).getData(sim.accumulators[1].ERROR.index);
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
    public static class VirialCO2TraPPEParam extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 353;
        public long numSteps = 1000;
        public double sigmaHSRef = 7;
        // the following??
       // composition
        public int[] nTypes = new int[]{1,1};
        public double refFrac = -1;
    }
}
