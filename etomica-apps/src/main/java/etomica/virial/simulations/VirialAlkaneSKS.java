/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.P2LennardJones;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
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

import java.awt.*;

/**
 * Mayer sampling simulation for alkanes using the SKS force field.
 */
public class VirialAlkaneSKS {


    public static void main(String[] args) {
        VirialSiepmannSpheresParam params = new VirialSiepmannSpheresParam();
        if (args.length > 0) {
            ReadParameters paramReader = new ReadParameters(args[0], params);
            paramReader.readParameters();
        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFrac = params.refFrac;
        double sigmaCH2 = 3.93;
        double sigmaCH3 = 3.93;
        double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;
        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        PotentialGroup pTargetGroup = new PotentialGroup(2);
        System.out.println("SKS "+nSpheres+"-mer chains B"+nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        double epsilonCH2 = Kelvin.UNIT.toSim(47.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(114.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        SpeciesAlkane species = new SpeciesAlkane(space, nSpheres);
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,
                          temperature,refCluster,targetCluster, nSpheres > 2);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);

        AtomType typeCH3 = species.getAtomType(0);
        AtomType typeCH2 = species.getAtomType(1);
        pTargetGroup.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));
        // CH2 on molecule1 to CH3 on molecule2
        pTargetGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));
        pTargetGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH2}));
        pTargetGroup.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));
        
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
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(0)});
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
            double size = (nSpheres+5)*1.5;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
            
            
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
            diameterManager.setDiameter(typeCH2, 0.8*sigmaCH2);
            diameterManager.setDiameter(typeCH3, 0.8*sigmaCH3);
            displayBox1.setDiameterHash(diameterManager);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
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
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }
        
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        
        sim.setAccumulatorBlockSize((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));
        if (nSpheres > 3) {
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
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(HSB[nPoints]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialSiepmannSpheresParam extends ParameterBase {
        public int nPoints = 2;
        public int nSpheres = 6;
        public double temperature = 300.0;   // Kelvin
        public long numSteps = 10000;
        public double refFrac = -1;
    }
}
