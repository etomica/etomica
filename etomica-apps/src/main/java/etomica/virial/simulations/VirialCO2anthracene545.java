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
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

/**
 * virial coefficients calculation
 * mixture of CO2 and anthracene
 * 
 * the code is modified from virialco2SKSmix
 * Mayer sampling to evaluate cluster integrals
 * CO2 is modeled as single site, anthracene is model as three site model 
 * this is from Iwai et al, Monte carlo simulation of solubilities of Na, PH and An supercritical, Fluid Phase Equilibria 1998
 * they modeled using 464 model or 545 model, this is 545 model
 * @author shu
 * Mar.4.2011
 * 
 */
public class VirialCO2anthracene545 {

    public static void main(String[] args) {
    	VirialCO2anthracene545Param params = new VirialCO2anthracene545Param();
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
        
        System.out.println("CO2(single) + anthracene(545 model) B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
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
        
        // cluster diagram of ref (hard sphere) system 
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        
        // CO2  potential, single site model, apply P2LennardJones
        double sigmaCO2 = 3.912;//they give in nm, need to be converted to angstrom
        double epsilonCO2 = Kelvin.UNIT.toSim(225.3);
        P2LennardJones pCO2 = new P2LennardJones(space, sigmaCO2, epsilonCO2);
        MayerGeneralSpherical fCO2 = new MayerGeneralSpherical(pCO2);
        MayerESpherical eCO2 = new MayerESpherical(pCO2);
    
  
               
        //this is  an-an potential 
        //there are two kinds of benzene rings. i.e. two kinds of interaction sites
        // name them as C and CH, C is "4"and CH is "5" 
        double sigmaC = 5.502;
        double sigmaCH = 5.502;
        double sigmaCCH = 5.502;
        double epsilonC = Kelvin.UNIT.toSim(137.3316);
        double epsilonCH = Kelvin.UNIT.toSim(250.9935);
        double epsilonCCH = Math.sqrt(epsilonC*epsilonCH);
        P2LennardJones p2C = new P2LennardJones(space, sigmaC, epsilonC);
        P2LennardJones p2CH = new P2LennardJones(space, sigmaCH, epsilonCH);
        P2LennardJones pCCH = new P2LennardJones(space,sigmaCCH , epsilonCCH);
        
       //make a potential group for anthracene
        
        PotentialGroup pAn = new PotentialGroup(2);
        MayerGeneral fAn= new MayerGeneral(pAn);
        MayerEGeneral eAn = new MayerEGeneral(pAn);
        
          // this is CO2-an potential, the 1st is CO2, 2nd is an
        double sigmaC_C = 4.707;
        double sigmaC_CH = 4.707;
        double epsilonC_C = Kelvin.UNIT.toSim(175.9);
        double epsilonC_CH = Kelvin.UNIT.toSim(237.8);

        P2LennardJones pC_C = new P2LennardJones(space, sigmaC_C, epsilonC_C);
        P2LennardJones pC_CH = new P2LennardJones(space,sigmaC_CH, epsilonC_CH);
        // make a potential group for co2 and anthracene
        PotentialGroup pCO2An = new PotentialGroup(2);
        MayerGeneral fCO2An= new MayerGeneral(pCO2An);
        MayerEGeneral eCO2An = new MayerEGeneral(pCO2An);
        // define target cluster, here we use mixture's cluster
        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fCO2,fCO2An},{fCO2An,fAn}},
                new MayerFunction[][]{{eCO2,eCO2An},{eCO2An,eAn}}, nTypes);
        targetCluster.setTemperature(temperature);
        

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		// initialize
      
      
        // put SpeciesFactoryCO2 and SpeciesFactoryAn here and these classes contain conformations already
        
        // this is CO2 species factory, apply spherical species factory, SpeciesFactorySpheres is a subclass of SpeciesFactory
        SpeciesFactorySpheres factoryCO2 = new SpeciesFactorySpheres();
        
        // this is anthracene species factory, apply  species factory
        SpeciesFactory factoryAn = new SpeciesFactory() {
            public ISpecies makeSpecies(Space space) {
            	SpeciesAnthracene3site545 species = new SpeciesAnthracene3site545(space);
                      return species;
            }
        };
              
        
        // now let us do the simulation!!!
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactory[]{factoryCO2,factoryAn}, nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},false);
        
        //put the species in the box
        SpeciesSpheresMono speciesCO2 = (SpeciesSpheresMono)sim.getSpecies(0);
        SpeciesAnthracene3site545 speciesAn = (SpeciesAnthracene3site545)sim.getSpecies(1);
        sim.integratorOS.setNumSubSteps(1000);

        AtomType typeC = speciesAn.getCType();
        AtomType typeCH = speciesAn.getCHType();
        AtomType typeCO2 = speciesCO2.getLeafType();

       // interaction between one site potential to the another site potential
        //between two solutes
        pAn.addPotential(p2C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));
        pAn.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH}));
        pAn.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeCH}));
        pAn.addPotential(pCCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH, typeC}));//switch
        
        //between solute and solvent
        pCO2An.addPotential(pC_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCO2, typeC}));
        pCO2An.addPotential(pC_CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCO2, typeCH}));

        //graphic part
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
    public static class VirialCO2anthracene545Param extends ParameterBase {
        public int nPoints = 5;
        public double temperature = 328.15;
        public long numSteps = 10000;
        public double sigmaHSRef = 7;
       //composition
        public int[] nTypes = new int[]{2,3};
        public double refFrac = -1;
    }
}
