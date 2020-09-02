/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.AlkaneEH.SpeciesMethane;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
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
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.cluster.VirialDiagramsMix2;

import java.awt.*;

/**
 *   Mayer sampling simulation for CO2(rigid, TraPPE)-CH4(rigid, TraPPE-EH) mixture
 *   cross virial coefficients
 *   Using VirialDiagramMix2 to generate diagrams
 *   rigid diagrams only since both of the components are rigid
 *   
 *   @author shu
 *   May 2013
 * 
 */
public class VirialCO2CH4Mix {
	
	public static void main(String[] args) {
        VirialMixParam params = new VirialMixParam();
        if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
			
		}
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int[] nTypes = params.nTypes;// composition
        double refFrac = params.refFrac;
        
        if ( nTypes[0]==0 || nTypes[1]==0){
        	throw new RuntimeException("refer to pure component virial coefficient calculation!");
        }
        if ( (nTypes[0]+nTypes[1])!= nPoints ){
        	throw new RuntimeException("wrong composition!");
        }
        System.out.println("\nCO2(TraPPE)+ CH4(TraPPE-EH) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
        System.out.println("both components are rigid, rigid diagram only");
        temperature = Kelvin.UNIT.toSim(temperature);
        boolean[] flex = new boolean[]{false,false};// both are rigid
                
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
        // ------------ ref cluster ------------------------- //
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
        //-------------- CO2 potential & CO2 mayer function-------------//
        P2CO2TraPPE pCO2 = new P2CO2TraPPE(space);
        MayerGeneral fCO2 = new MayerGeneral(pCO2);
        // ------------- CH4 potential & CH4 mayer function---------- //
        double sigmaH = 3.31;//middle point of C-H bond
        double sigmaC = 3.31;//carbon core
        double sigmaCH = 3.31;
        double epsilonH = Kelvin.UNIT.toSim(15.3);
        double epsilonC = Kelvin.UNIT.toSim(0.01);
        double epsilonCH = Math.sqrt((epsilonH * epsilonC));
        P2LennardJones p2HH = new P2LennardJones(space, sigmaH, epsilonH);
        P2LennardJones p2CC = new P2LennardJones(space, sigmaC, epsilonC);
        P2LennardJones p2CH = new P2LennardJones(space, sigmaCH, epsilonCH);
        PotentialGroup pCH4 = new PotentialGroup(2);
        MayerGeneral fCH4= new MayerGeneral(pCH4);
        // ------------ CO2(1st)-CH4(2nd) potential & mayer function ------------//
        P2LennardJones pC_H = new P2LennardJones(space, 0.5*(sigmaH + pCO2.getSigmaC()), Math.sqrt(epsilonH * pCO2.getEpsilonC()));
        P2LennardJones pC_C = new P2LennardJones(space, 0.5*(sigmaC + pCO2.getSigmaC()), Math.sqrt(epsilonC * pCO2.getEpsilonC()));
        P2LennardJones pO_H = new P2LennardJones(space, 0.5*(sigmaH + pCO2.getSigmaO()), Math.sqrt(epsilonH * pCO2.getEpsilonO()));
        P2LennardJones pO_C = new P2LennardJones(space, 0.5*(sigmaC + pCO2.getSigmaO()), Math.sqrt(epsilonC * pCO2.getEpsilonO()));
        PotentialGroup pCO2CH4 = new PotentialGroup(2);
        MayerGeneral fCO2CH4= new MayerGeneral(pCO2CH4); 
        
        // use VirialDiagram to construct target cluster
        VirialDiagramsMix2 diagrams = new VirialDiagramsMix2(nPoints,flex);
        diagrams.setDoReeHoover(false);
        MayerFunction[][] f = new MayerFunction[][]{{fCO2,fCO2CH4},{fCO2CH4,fCH4}};
        ClusterSum targetCluster = diagrams.makeVirialCluster(f,-1, nTypes);//flexID always -1 since both are rigid
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        targetDiagrams = diagrams.makeSingleVirialClusters(targetCluster, f, -1, nTypes);//flexID always -1 since both are rigid
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        SpeciesGeneral speciesCO2 = SpeciesTraPPECO2.create(space);
        SpeciesGeneral speciesCH4 = SpeciesMethane.create(false);

        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), 
        		                                             ClusterWeightAbs.makeWeightCluster(targetCluster)};

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesCO2,speciesCH4}, 
        		nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},targetDiagrams,sampleClusters,false);

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        AtomType typeC_CH4 = speciesCH4.getAtomType(0);//C in CH4
        AtomType typeH = speciesCH4.getAtomType(1);//H in CH4
        AtomType typeC_CO2 = speciesCO2.getAtomType(0);//  C in CO2
        AtomType typeO = speciesCO2.getAtomType(1);// O in CO2
        
        // CH4 potential
        pCH4.addPotential(p2HH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));
        pCH4.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CH4, typeH}));// H on molecule1 --- C on molecule2
        pCH4.addPotential(p2CH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CH4, typeH}));// C on molecule1 --- H on molecule2
        pCH4.addPotential(p2CC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CH4, typeC_CH4}));
        // CO2-CH4 potential
        pCO2CH4.addPotential(pC_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CO2, typeH}));
        pCO2CH4.addPotential(pC_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CO2, typeC_CH4}));
        pCO2CH4.addPotential(pO_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeH}));
        pCO2CH4.addPotential(pO_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeC_CH4}));
                
        sim.integratorOS.setNumSubSteps(1000);

        // find a proper configuration 
        double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        for (int j=0; j<10000 && (pi < 1e-10 || Double.isNaN(pi)); j++) {
            sim.integrators[1].doStep();
            pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        }
        if ( pi == 0) {
            throw new RuntimeException("could not find a configuration for target system");
        }
        sim.accumulators[1].reset();// don't want to collect these data!!!!
        
        if(false) {
    double size = 10;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox dBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox dBox1 = simGraphic.getDisplayBox(sim.box[1]);
            dBox0.setPixelUnit(new Pixel(300.0 / size));
            dBox1.setPixelUnit(new Pixel(300.0 / size));
            dBox0.setShowBoundary(false);
            dBox1.setShowBoundary(false);

            //set diameters
            DiameterHashByType diameter = new DiameterHashByType();
            diameter.setDiameter(speciesCO2.getAtomType(0), 0.2);
            diameter.setDiameter(speciesCO2.getAtomType(1), 0.3);
            diameter.setDiameter(speciesCH4.getAtomType(0), 0.3);
            diameter.setDiameter(speciesCH4.getAtomType(1), 0.4);

            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(speciesCO2.getAtomType(0), Color.blue);
            colorScheme.setColor(speciesCO2.getAtomType(1), Color.red);
            colorScheme.setColor(speciesCH4.getAtomType(0), Color.green);
            colorScheme.setColor(speciesCH4.getAtomType(1), Color.yellow);

            ((DisplayBoxCanvasG3DSys) dBox1.canvas).setBackgroundColor(Color.WHITE);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(100);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
    sim.equilibrate(null, 20, false);
    sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
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
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
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
                    if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }
sim.getController().runActivityBlocking(ai);

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
        
	}

    /**
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 5;
        public double temperature = 298;
        public long numSteps = 1000000;
        public double sigmaHSRef = 4.5;
        public int[] nTypes = new int[]{3,2};
        public double refFrac = -1;
    }
}
