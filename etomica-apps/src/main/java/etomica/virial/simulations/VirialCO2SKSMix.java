/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Oxygen;
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
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialCO2SKSMix {


    public static void main(String[] args) {
        VirialCO2Param params = new VirialCO2Param();
        if (args.length > 0) {
            ReadParameters readParameters = new ReadParameters(args[0], params);
            readParameters.readParameters();
        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int nSpheres = params.nSpheres;
        int[] nTypes = params.nTypes;
        double refFrac = params.refFrac;
        int sum = 0;
        for (int i=0; i<nTypes.length; i++) {
            if (nTypes[i] == 0) {
                throw new RuntimeException("for pure component, use a different class");
            }
            sum += nTypes[i];
        }
        if (sum != nPoints) {
            throw new RuntimeException("Number of each type needs to add up to nPoints");
        }
        
        System.out.println("CO2 + AlkaneSKS overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
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
        P2CO2EMP2 pCO2 = new P2CO2EMP2(space);
        MayerGeneral fCO2 = new MayerGeneral(pCO2);
        MayerEGeneral eCO2 = new MayerEGeneral(pCO2);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        
        double sigmaCH2 = 3.93;
        double sigmaCH3 = 3.93;
        double epsilonCH2 = Kelvin.UNIT.toSim(47.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(114.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        PotentialGroup pAlkane = new PotentialGroup(2);

        MayerGeneral fAlkane= new MayerGeneral(pAlkane);
        MayerEGeneral eAlkane = new MayerEGeneral(pAlkane);

        P2LennardJones pC_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pCO2.getSigmaC()), Math.sqrt(epsilonCH2*pCO2.getEpsilonC()));
        P2LennardJones pC_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pCO2.getSigmaC()), Math.sqrt(epsilonCH3*pCO2.getEpsilonC()));
        P2LennardJones pO_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pCO2.getSigmaO()), Math.sqrt(epsilonCH2*pCO2.getEpsilonO()));
        P2LennardJones pO_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pCO2.getSigmaO()), Math.sqrt(epsilonCH3*pCO2.getEpsilonO()));
        
        PotentialGroup pCO2Alkane = new PotentialGroup(2);
        MayerGeneral fCO2Alkane= new MayerGeneral(pCO2Alkane);
        MayerEGeneral eCO2Alkane = new MayerEGeneral(pCO2Alkane);

        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fCO2,fCO2Alkane},{fCO2Alkane,fAlkane}},
                new MayerFunction[][]{{eCO2,eCO2Alkane},{eCO2Alkane,eAlkane}}, nTypes);
        targetCluster.setTemperature(temperature);
        

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		
        final IConformation conformation = new IConformation() {
            
            public void initializePositions(IAtomList atomList) {
                // atoms are C, O and O, so we arrange them as 1-0-2
                double bondL = 1.1491;
                atomList.get(0).getPosition().E(0);
                atomList.get(1).getPosition().E(0);
                atomList.get(1).getPosition().setX(0, -bondL);
                atomList.get(2).getPosition().E(0);
                atomList.get(2).getPosition().setX(0, +bondL);
            }
        };
        SpeciesSpheresHetero speciesCO2 = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
        speciesCO2.setChildCount(new int[]{1,2});
        speciesCO2.setConformation(conformation);

        SpeciesAlkane speciesAlkane = new SpeciesAlkane(space, nSpheres);
        
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesCO2,speciesAlkane}, nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},true);
        ((MCMoveClusterWiggleMulti)sim.mcMoveWiggle[0]).setSpecies(sim.getSpecies(1));
        ((MCMoveClusterWiggleMulti)sim.mcMoveWiggle[1]).setSpecies(sim.getSpecies(1));
        sim.integratorOS.setNumSubSteps(1000);

        AtomType typeCH3 = speciesAlkane.getAtomType(0);
        AtomType typeCH2 = speciesAlkane.getAtomType(1);
        AtomType typeC = speciesCO2.getAtomType(0);
        AtomType typeO = speciesCO2.getAtomType(1);
        pAlkane.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));
        // CH2 on molecule1 to CH3 on molecule2
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH2}));
        pAlkane.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));

        pCO2Alkane.addPotential(pC_CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH3}));
        pCO2Alkane.addPotential(pC_CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH2}));
        pCO2Alkane.addPotential(pO_CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeCH3}));
        pCO2Alkane.addPotential(pO_CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeCH2}));

        
        sim.integratorOS.setNumSubSteps(1000);

        // create the intramolecular potential here, add to it and add it to
        // the potential master if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        if (nSpheres > 2) {
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
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres-3][4];
            for (int i=0; i<nSpheres-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            pIntra.addPotential(p4, new Atomset4IteratorIndexList(quads));
            torsionMoves = new MCMoveClusterTorsionMulti[2];
            torsionMoves[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves[0].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMoves[0]);
            torsionMoves[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves[1].setTemperature(temperature);
            sim.integrators[1].getMoveManager().addMCMove(torsionMoves[1]);
        }
        if (nSpheres > 4) {
            pIntra.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres-1}}));
        }
        if (nSpheres > 5) {
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

        if (false) {
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
            DiameterHashByType diameterHash = (DiameterHashByType) dBox1.getDiameterHash();
            diameterHash.setDiameter(((SpeciesAlkane) sim.getSpecies(1)).getCH2Type(), sigmaCH2);
            diameterHash.setDiameter(((SpeciesAlkane) sim.getSpecies(1)).getCH3Type(), sigmaCH3);
            dBox0.setDiameterHash(diameterHash);
            ColorSchemeByType colorScheme = (ColorSchemeByType) dBox1.getColorScheme();
            colorScheme.setColor(((SpeciesAlkane) sim.getSpecies(1)).getCH2Type(), new Color(190, 190, 190));
            colorScheme.setColor(((SpeciesAlkane) sim.getSpecies(1)).getCH3Type(), new Color(240, 240, 240));
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), new Color(100, 100, 150));
            colorScheme.setColor(sim.getSpecies(0).getAtomType(1), Color.RED);
            dBox0.setColorScheme(colorScheme);
            ((DisplayBoxCanvasG3DSys) dBox1.canvas).setBackgroundColor(Color.WHITE);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController2().addActivity(new ActivityIntegrate2(sim.integratorOS));
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
    public static class VirialCO2Param extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 300;
        public long numSteps = 100000;
        public double sigmaHSRef = 7;
        public int nSpheres = 6;
        public int[] nTypes = new int[]{1,1};
        public double refFrac = -1;
    }
}
