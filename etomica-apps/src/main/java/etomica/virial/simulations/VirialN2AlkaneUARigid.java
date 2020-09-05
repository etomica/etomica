/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Nitrogen;
import etomica.config.IConformation;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
  * cross virial coefficients of mixture of N2 and alkanes(TraPPE-UA)
  * no flexible diagrams, just rigid diagrams
  * just for B11, see "VirialN2AlkaneMix" for flexible corrections
  * N2:use CO2EMP, C<==> "A", O<==>"N"
  * 
  * @author shu
  * Feb-2013
 */
public class VirialN2AlkaneUARigid {


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
        double sigmaCH2= params.sigmaCH2;
        double sigmaCH3= params.sigmaCH3;
        int nSpheres = params.nSpheres;
       
        int[] nTypes = params.nTypes;
        double refFrac = params.refFrac;
        if (nSpheres ==1){
        	throw new RuntimeException("I cannot handle methane!");
        }
        if ( (nTypes[0] != 1) || (nTypes[1]!=1) ){
        	throw new RuntimeException("I can only calculate B11!");
        }
        System.out.println("rigid diagrams only");
        System.out.println("N2(TraPPE)+"+nSpheres+"smer n-alkane(TraPPE-UA) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
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
        // ------------ potential parameters ------------------------- //
        double sigmaA = 1.0;// set arbitrarily
        double sigmaN = 3.31;
        double epsilonA = 0.0;
        double epsilonN = 36.0;
        double charge = 0.964;// in the center of mass
        double sigmaAN= (sigmaA + sigmaN) * 0.5 ;
        double epsilonAN = 0.0;
        //-------------- N2 potential -------------//
        P2CO2EMP pN2 = new P2CO2EMP(space, sigmaA, sigmaAN, sigmaN,epsilonA ,epsilonAN, Kelvin.UNIT.toSim(epsilonN), Electron.UNIT.toSim(charge));
        MayerGeneral fN2 = new MayerGeneral(pN2);
        MayerEGeneral eN2 = new MayerEGeneral(pN2);
        // ------------- Alkane potential ---------- //
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        PotentialGroup pAlkane = new PotentialGroup(2);

        MayerGeneral fAlkane= new MayerGeneral(pAlkane);
        MayerEGeneral eAlkane = new MayerEGeneral(pAlkane);
        // N2-alkane potential, the 1st is N2, the 2nd is alkane
 //       P2LennardJones pA_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pN2.getSigmaC()), Math.sqrt(epsilonCH2*pN2.getEpsilonC()));
 //       P2LennardJones pA_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pN2.getSigmaC()), Math.sqrt(epsilonCH3*pN2.getEpsilonC()));
        P2LennardJones pN_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pN2.getSigmaO()), Math.sqrt(epsilonCH2*pN2.getEpsilonO()));
        P2LennardJones pN_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pN2.getSigmaO()), Math.sqrt(epsilonCH3*pN2.getEpsilonO()));
        PotentialGroup pN2Alkane = new PotentialGroup(2);
        MayerGeneral fN2Alkane= new MayerGeneral(pN2Alkane);
        MayerEGeneral eN2Alkane = new MayerEGeneral(pN2Alkane);

        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fN2,fN2Alkane},{fN2Alkane,fAlkane}},
                new MayerFunction[][]{{eN2,eN2Alkane},{eN2Alkane,eAlkane}}, nTypes);
        targetCluster.setTemperature(temperature);

        AtomType nType = AtomType.element(Nitrogen.INSTANCE);
        double bondL = 1.10;
        SpeciesGeneral speciesN2 = new SpeciesBuilder(space)
                .addAtom(AtomType.simple("A"), space.makeVector())
                .addAtom(nType, Vector.of(-0.5 * bondL))
                .addAtom(nType, Vector.of(+0.5 * bondL))
                .build();
        SpeciesGeneral speciesAlkane = SpeciesAlkane.create(nSpheres);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesN2,speciesAlkane}, nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},false);
       
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;
        
        MCMoveClusterWiggleMulti[] wiggleMove = new MCMoveClusterWiggleMulti[2];

        // if ethane, then do not do wiggle
        if (nSpheres==2){
        	
        	System.out.println("n=2, no wiggles !");

        }
        else {
        	System.out.println("not ethane, add wiggle move to alkane !");
        	for (int i=0; i<2; i++) {
        		wiggleMove[i] = new MCMoveClusterWiggleMulti(sim, sim.integrators[i].getPotentialMaster(), targetCluster.pointCount(), space);
        		wiggleMove[i].setSpecies(sim.getSpecies(1));
        		sim.integrators[i].getMoveManager().addMCMove(wiggleMove[i]);
        	}
       
        }
        
        sim.integratorOS.setNumSubSteps(1000);
        // alkane potential
        AtomType typeCH3 = speciesAlkane.getAtomType(0);
        AtomType typeCH2 = speciesAlkane.getAtomType(1);
        AtomType typeA = speciesN2.getAtomType(0);
        AtomType typeN = speciesN2.getAtomType(1);
        pAlkane.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));
        // CH2 on molecule1 to CH3 on molecule2
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH2}));
        pAlkane.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));
        // N2-alkane potential
        //pN2Alkane.addPotential(pA_CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeA, typeCH3}));
        //pN2Alkane.addPotential(pA_CH2, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeA, typeCH2}));
        pN2Alkane.addPotential(pN_CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeCH3}));
        pN2Alkane.addPotential(pN_CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeCH2}));
        
        sim.integratorOS.setNumSubSteps(1000);

        // create the intramolecular potential here, add to it and add it to the potential master if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        //*********************************************************************************************************//
        // intramolecular potential
        //*********************************************************************************************************//
        if (nSpheres > 2) {
        	System.out.println("n>2, add bending potential to intramolcular potential!");
            P3BondAngle p3 = new P3BondAngle(space);
            p3.setAngle(Math.PI*114.0/180.0);
            p3.setEpsilon(Kelvin.UNIT.toSim(62500));
            int[][] triplets = new int[nSpheres-2][3];
            for (int i=0; i<nSpheres-2; i++) {
                triplets[i][0] = i;
                triplets[i][1] = i+1;
                triplets[i][2] = i+2;
            }
            pIntra.addPotential(p3, new Atomset3IteratorIndexList(triplets));
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(1)});
        }
        MCMoveClusterTorsionMulti[] torsionMoves = null;
        if (nSpheres > 3) {
        	System.out.println("n>3, add torsion potential to intramolecular potential!");
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres-3][4];
            for (int i=0; i<nSpheres-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            
            pIntra.addPotential(p4, new Atomset4IteratorIndexList(quads));
            
            System.out.println("n>3, add torsion mc move!");
            torsionMoves = new MCMoveClusterTorsionMulti[2];
            for (int j = 0; j<2; j++){   
            	torsionMoves[j] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            	torsionMoves[j].setTemperature(temperature);
                torsionMoves[j].setSpecies(sim.getSpecies(1));
                sim.integrators[j].getMoveManager().addMCMove(torsionMoves[j]);

            }
            
        }
        if (nSpheres > 4) {
        	System.out.println("n>4, add CH3-CH3 pair intramolecular potential!");
            pIntra.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres-1}}));
        }
        if (nSpheres > 5) {
        	System.out.println("n>5, add CH3-CH2 pair intramolecular potential!");
            int[][] pairs = new int[2*(nSpheres-5)][2];
            for (int i=0; i<nSpheres-5; i++) {
                pairs[2*i][0] = 0;
                pairs[2*i][1] = nSpheres-2-i;
                pairs[2*i+1][0] = nSpheres-1;
                pairs[2*i+1][1] = i+1;
            }
            pIntra.addPotential(p2CH2CH3,new ApiIndexList(pairs));
        }
        if (nSpheres > 6) {
        	System.out.println("n>6, add CH2-CH2 pair intramolecular potential!");
            int[][] pairs = new int[(nSpheres-6)*(nSpheres-5)/2][2];
            int k = 0;
            for (int i=1; i<nSpheres-5; i++) {
                for (int j=i+4; j<nSpheres-1; j++) {
                    pairs[k][0] = i;
                    pairs[k][1] = j;
                    k++;
                }
            }
            pIntra.addPotential(p2CH2,new ApiIndexList(pairs));
        }

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
            diameter.setDiameter(speciesN2.getAtomType(0), 0.2);
            diameter.setDiameter(speciesN2.getAtomType(1), 0.3);
            diameter.setDiameter(speciesAlkane.getTypeByName("CH2"), 0.3);
            diameter.setDiameter(speciesAlkane.getTypeByName("CH3"), 0.4);

            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(speciesN2.getAtomType(0), Color.blue);
            colorScheme.setColor(speciesN2.getAtomType(1), Color.red);
            colorScheme.setColor(speciesAlkane.getTypeByName("CH2"), Color.green);
            colorScheme.setColor(speciesAlkane.getTypeByName("CH3"), Color.yellow);

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
        if(nSpheres>2){
        	//print step size of wiggle
        	System.out.println("Wiggle move step sizes " + wiggleMove[0].getStepSize() + " "+	wiggleMove[1].getStepSize());
        }
        if (nSpheres > 3) { // print torsion ratio
            System.out.println("Torsion move acceptance "+torsionMoves[0].getTracker().acceptanceRatio()+" "+
                    torsionMoves[1].getTracker().acceptanceRatio());
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
        public int nPoints = 2;
        public double temperature = 298;
        public long numSteps = 1000000;
        public int nSpheres = 8;
        public double sigmaCH2 = 3.95;// in alkane
        public double sigmaCH3 = 3.75;// in alkane
        public double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;

        public int[] nTypes = new int[]{1,1};
        public double refFrac = -1;
    }
}
