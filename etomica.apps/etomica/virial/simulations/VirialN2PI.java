/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import java.util.List;

import etomica.api.IAtom;
import etomica.api.IAtomType;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IMoleculeList;
import etomica.atom.AtomHydrogen;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.IAtomTypeOriented;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Nitrogen;
import etomica.config.ConformationLinear;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.P2NitrogenHellmann;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.ClusterWheatleySoft;
import etomica.virial.MCMoveClusterRingRegrow;
import etomica.virial.MCMoveClusterRingRegrowOrientation;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.PotentialGroupPI;
import etomica.virial.cluster.Standard;

public class VirialN2PI {
    public static void main(String[] args) {
        VirialN2Param params = new VirialN2Param();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        else {
            // default options - choose these before committing to CVS
            params.nBeads = 8;
            params.temperatureK = 500;
            params.numSteps = (long)1E6;
            params.pN2HellmannA = false;

            // runtime options - make changes in these and not the default options above           
//            params.potentialLevel = level.semiClassical;
//            params.temperatureK = 500;
//            params.numSteps = (long)1E8;
        }
        final int nPoints = params.nPoints;
        final double temperatureK = params.temperatureK;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final boolean pairOnly = params.pairOnly;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        final boolean pN2HellmannA = params.pN2HellmannA;
        final int nBeads = params.nBeads;
        final int beadFac = params.beadFac;        
        

        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        
        System.out.println("Overlap sampling for N2 pair potential of Hellmann (2013) at " + temperatureK + " K");
        System.out.println("Path Integral Monte Carlo (PIMC) calculation with P = "+nBeads);

        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
        
        Space space = Space3D.getInstance();
        // make ref and tar clusters
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);
        
        final P2NitrogenHellmann p2 = new P2NitrogenHellmann(space);
        if (pN2HellmannA) p2.parametersB = false;
        PotentialGroupPI pTarGroup = new PotentialGroupPI(beadFac);
        pTarGroup.addPotential(p2, new ApiIntergroupCoupled());
        MayerGeneral fTar = new MayerGeneral(pTarGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {                
                return super.f(pair, r2, beta/nBeads);
            }
        };
        ClusterWheatleySoft tarCluster = new ClusterWheatleySoft(nPoints, fTar, 1e-12);
        tarCluster.setTemperature(temperature);
        
        // make species        
        IAtomTypeOriented atype = new AtomTypeOrientedSphere(Nitrogen.INSTANCE, space);
        SpeciesSpheresHetero speciesN2 = null;        
        speciesN2 = new SpeciesSpheresHetero(space,new IAtomTypeOriented [] {atype}) {
            protected IAtom makeLeafAtom(IAtomType leafType) {
                return new AtomHydrogen(space,(IAtomTypeOriented)leafType,blN2);
            }
        };
        speciesN2.setChildCount(new int [] {nBeads});
        speciesN2.setConformation(new ConformationLinear(space, 0));        
        // make simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesN2, temperature, refCluster, tarCluster);
//        sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;
        
//        add additional moves here, simulation already has translation and rotation moves
//        rotation is a bit pointless when we can regrow the ring completely

        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
            sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        }        
        
//        ring regrow translation
        MCMoveClusterRingRegrow refTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*massN*temperature);
        refTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

        MCMoveClusterRingRegrow tarTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
        tarTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().addMCMove(refTr);
            sim.integrators[1].getMoveManager().addMCMove(tarTr);
        }
//        ring regrow orientation
        MCMoveClusterRingRegrowOrientation refOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);        
        MCMoveClusterRingRegrowOrientation tarOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
        refOr.setStiffness(temperature, massN);
        tarOr.setStiffness(temperature, massN);
        
        if (nBeads != 1) {
            sim.integrators[0].getMoveManager().addMCMove(refOr);
            sim.integrators[1].getMoveManager().addMCMove(tarOr);
        }
        
        System.out.println();
        String refFileName = null;
        if (isCommandline) {
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints;
            refFileName += pairOnly ? "_2b" : "_3b";
            refFileName += "_"+tempString+"_PI";
        }
        long t1 = System.currentTimeMillis();
        
        sim.initRefPref(refFileName, steps/40);
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }
        sim.equilibrate(refFileName, steps/10);
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.ai.setMaxSteps(1000);
        
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());
        
        final double refIntegralF = HSB[nPoints];
        if (! isCommandline) {
            IIntegratorListener progressReport = new IIntegratorListener() {        
                public void integratorInitialized(IIntegratorEvent e) {}
                public void integratorStepStarted(IIntegratorEvent e) {}
                public void integratorStepFinished(IIntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    //                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }
        // this is where the simulation takes place
        sim.getController().actionPerformed();
        //end of simulation
        long t2 = System.currentTimeMillis();
        
        
        System.out.println(" ");
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        System.out.println(" ");
        
        System.out.println("Reference system: ");
        List<MCMove> refMoves = sim.integrators[0].getMoveManager().getMCMoves();
        for (MCMove m : refMoves) {
            double acc = m.getTracker().acceptanceRatio();
            System.out.println(m.toString()+" acceptance ratio: "+acc);            
        }
        System.out.println("Target system: ");
        List<MCMove> tarMoves = sim.integrators[1].getMoveManager().getMCMoves();
        for (MCMove m : tarMoves) {
            double acc = m.getTracker().acceptanceRatio();

//            if (acc == 1) {
//                throw new RuntimeException("something seems fishy");
//            }
            System.out.println(m.toString()+" acceptance ratio: "+acc);
        }
        sim.printResults(HSB[nPoints]);        
        
        if ((t2-t1)/1000.0 > 24*3600) {
            System.out.println("time: "+(t2-t1)/(24*3600*1000.0)+" days");
        }
        else if ((t2-t1)/1000.0 > 3600) {
            System.out.println("time: "+(t2-t1)/(3600*1000.0)+" hrs");
        }
        else if ((t2-t1)/1000.0 > 60) {
            System.out.println("time: "+(t2-t1)/(60*1000.0)+" mins");
        }
        else {
            System.out.println("time: "+(t2-t1)/1000.0+" secs");
        }
    }    
    public static final double massN = Nitrogen.INSTANCE.getMass();
    public static final double blN2 = 1.1014;
    /**
     * Inner class for parameters
     */
    public static class VirialN2Param extends ParameterBase {
        public int nPoints = 2;
        public int nBeads = 8;
        public double temperatureK = 500.0;   // Kelvin
        public long numSteps = 1000000;
        public double refFrac = -1;        
        public double sigmaHSRef = 4.50; // -1 means use equation for sigmaHSRef
        public boolean pairOnly = true;
        public boolean pN2HellmannA = true;
        public int beadFac = 2;        
    }

}
