/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterRadiusGyration;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeRandom;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.compute.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.mcmove.MCMoveClusterAngle;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.io.FileWriter;
import etomica.zeno.MeterZENO;
import etomica.zeno.WalkerExterior;

/**
 * Single-molecule calculations for input into theta temperature calculations
 */
public class VirialChainSingle {

    private static int[] MeterRadiusGyration;

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nSpheres = 128;
            params.temperature = 4; //.42304894079030;
            params.kBend = 0;
            params.numSteps = 1000000;
        }
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double eFENE = params.eFENE;
        double rc = params.rc;
        double bondLength = params.bondLength;
        double kBend = params.kBend;
        boolean ideal = params.ideal;

        if (ideal) {
            eFENE = 0;
            temperature = 1;
            kBend = 0;
        }

        Space space = Space3D.getInstance();

        ConformationLinear conf = new ConformationLinear(Space3D.getInstance(), bondLength, new double[]{0,Math.PI/2});
        AtomType type = AtomType.simple("A");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(type, nSpheres)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        if (ideal) {
            pTarget  = new PotentialMoleculePair(space, sm) {
                public double energy(IMoleculeList molecules) {
                    return 0;
                }
            };
        }
        System.out.println(nSpheres+"-mer chain at T = "+temperature);
        System.out.println("Bond length: "+bondLength+"  with eFENE "+eFENE);
        if (kBend == Double.POSITIVE_INFINITY) System.out.println("Rigid bond angles");
        else if (kBend == 0) System.out.println("Flexible bond angles");
        else System.out.println("Bond angle bend constant: "+kBend);
        if (rc>0 && rc<Double.POSITIVE_INFINITY) System.out.println("LJ truncated at "+rc+" with "+params.truncation);
        else System.out.println("LJ untruncated");
        IPotential2 p2 = new P2LennardJones(1, 1);
        if (rc > 0 && rc < Double.POSITIVE_INFINITY) {
            if (params.truncation == TruncationChoice.SIMPLE) {
                p2 = new P2SoftSphericalTruncated(p2, rc);
            }
            else if (params.truncation == TruncationChoice.SHIFT){
                p2 = new P2SoftSphericalTruncatedShifted(p2, rc);
            }
            else if (params.truncation == TruncationChoice.FORCESHIFT) {
                p2 = new P2SoftSphericalTruncatedForceShifted(p2, rc);
            }
        }
        if (!ideal) pTarget.setAtomPotential(type, type, p2);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        if (nSpheres > 1 && (kBend == 0 || eFENE > 0)) {
            IPotential2 pBonding;
            if (eFENE > 0 && eFENE < Double.POSITIVE_INFINITY) {
                pBonding = new P2Fene(2, eFENE);
            }
            else {
                // we need to do this to convince the system that the molecules are not rigid
                // if bondingInfo thinks molecules are rigid then intramolecular LJ will not be computed
                pBonding = new IPotential2() {
                    @Override
                    public double getRange() { return 2; }
                    @Override
                    public void u012add(double r2, double[] u012) { }
                };
            }
            List<int[]> pairs = new ArrayList<>();
            for (int i=0; i<nSpheres-1; i++) {
                pairs.add(new int[]{i,i+1});
            }
            bondingInfo.setBondingPotentialPair(species, pBonding, pairs);
        }
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }
        if (kBend < Double.POSITIVE_INFINITY) {
            P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(kBend);
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }

        Simulation sim = new Simulation(Space3D.getInstance(), sm);
        sim.makeBox(new BoundaryRectangularNonperiodic(space));
        sim.box().setNMolecules(species, 1);

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box(), bondingInfodk);
        P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(1);
        bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);

        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, sim.box(), bondingInfo);
        NeighborManager nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo);
        if (ideal) {
            nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo) {
                public NeighborIterator makeNeighborIterator() {
                    return new NeighborIterator() {
                        @Override
                        public void iterUpNeighbors(int i, NeighborConsumer consumer) {}

                        @Override
                        public void iterDownNeighbors(int i, NeighborConsumer consumer) {}

                        @Override
                        public void iterAllNeighbors(int i, NeighborConsumer consumer) {}

                        @Override
                        public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {return 0;}
                    };
                }

            };
        }
        PotentialComputePair pcPair = new PotentialComputePair(sm, sim.box(), nbrManager, pTarget.getAtomPotentials());
        PotentialComputeAggregate pc = new PotentialComputeAggregate(pmBonding, pcPair);

        IntegratorMC integrator = new IntegratorMC(pc, sim.getRandom(), temperature, sim.box());

        System.out.println(steps+" steps");

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        MCMoveClusterAngle angleMove1 = null;
        MCMoveClusterReptate reptateMove = null;
        MCMoveClusterShuffle shuffleMove = null;
        if (kBend < Double.POSITIVE_INFINITY) {
            angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
            angleMove1.setBox(sim.box());
            integrator.getMoveManager().addMCMove(angleMove1);

            reptateMove = new MCMoveClusterReptate(pc, space, sim.getRandom());
            reptateMove.setBox(sim.box());
            integrator.getMoveManager().addMCMove(reptateMove);

            shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
            shuffleMove.setBox(sim.box());
            integrator.getMoveManager().addMCMove(shuffleMove);
            ((MCMoveStepTracker) shuffleMove.getTracker()).setAcceptanceTarget(0.3);
        }

        if (false) {
            ActivityIntegrate ai = new ActivityIntegrate(integrator);
            sim.getController().addActivity(ai, Long.MAX_VALUE, 10);

            double size = (nSpheres + 5) * 1.5;
            sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "stuff", 1);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            simGraphic.makeAndDisplayFrame();

            return;
        }

        long t1 = System.nanoTime();
        // if running interactively, don't use the file
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        System.out.println("equilibration finished");
        System.out.println("Angle move step size    " + angleMove1.getStepSize());
        System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);

        DataSourceThetaSingle dsTheta = new DataSourceThetaSingle(sim.box(), temperature, kBend, pcPair, pmBondingdk);
        AccumulatorAverageFixed accTheta = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpTheta = new DataPumpListener(dsTheta, accTheta, 10);
        integrator.getEventManager().addListener(pumpTheta);

        MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());

        AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        integrator.getEventManager().addListener(pumpRg);
        // ---- ZENO: Rh, intrinsic viscosity, charge distribution ----
        double[] zenoSigma = new double[1];
        zenoSigma[type.getIndex()] = 1.0;                  // the only atom type, diameter 1
        MeterZENO meterZeno = new MeterZENO(sim.box(), sim.getRandom(), zenoSigma, false, WalkerExterior.Shape.SPHERE);
        meterZeno.setEigenvaluesForPade(true);
        meterZeno.setNumWalks(params.zenoWalks);
        meterZeno.setChargeBinWidth(params.chargeBinWidth);
        AccumulatorAverageFixed accZeno = new AccumulatorAverageFixed(10);
        DataPumpListener pumpZeno = new DataPumpListener(meterZeno, accZeno, params.zenoInterval);
        integrator.getEventManager().addListener(pumpZeno);


        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();
        if (kBend < Double.POSITIVE_INFINITY) {
            System.out.println("Reptate move acceptance " + reptateMove.getTracker().acceptanceProbability());
            System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability());
            System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
            System.out.println();
        }

        DataGroup allYourBase = (DataGroup)accTheta.getData();
        IData avg = allYourBase.getData(accTheta.AVERAGE.index);
        IData err = allYourBase.getData(accTheta.ERROR.index);
        IData cor = allYourBase.getData(accTheta.BLOCK_CORRELATION.index);

        System.out.println("dlnqdb: "+avg.getValue(0)+"   err: "+err.getValue(0)+"   cor: "+cor.getValue(0));
        System.out.println("dlnqdk: "+avg.getValue(1)+"   err: "+err.getValue(1)+"   cor: "+cor.getValue(1));
        System.out.println("qkboq: "+avg.getValue(2)+"   err: "+err.getValue(2)+"   cor: "+cor.getValue(2));
        System.out.println("qbboq: "+avg.getValue(3)+"   err: "+err.getValue(3)+"   cor: "+cor.getValue(3));

        double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);

        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);


        // ---- ZENO results ----
        IData zenoAvg = accZeno.getData(accZeno.AVERAGE);
        IData zenoErr = accZeno.getData(accZeno.ERROR);
        IData zenoCor = accZeno.getData(accZeno.BLOCK_CORRELATION);
        double avgRh = zenoAvg.getValue(0), errRh = zenoErr.getValue(0), corRh = zenoCor.getValue(0);
        double avgEta = zenoAvg.getValue(5), errEta = zenoErr.getValue(5), corEta = zenoCor.getValue(5);
        System.out.println("Rh: " + avgRh + "   err: " + errRh + "   cor: " + corRh);
        System.out.println("eta: " + avgEta + "   err: " + errEta + "   cor: " + corEta);

        // results_chain.csv: one row per run
        String outFile = System.getProperty("user.dir") + "/results_chain.csv";
        boolean outExists = new java.io.File(outFile).exists();
        try (FileWriter fw = new FileWriter(outFile, true)) {
            if (!outExists) fw.write("nSpheres,temperature,ideal,rc,Rg2,Rg2_err,Rg2_cor,Rh,Rh_err,Rh_cor,eta,eta_err,eta_cor,numSteps,zenoWalks,zenoInterval\n");
            fw.write(String.format("%d,%.4f,%b,%.3f,%.8f,%.8f,%.6f,%.8f,%.8f,%.6f,%.8f,%.8f,%.6f,%d,%d,%d%n",
                    nSpheres, temperature, ideal, rc, avgRg, errRg, corRg,
                    avgRh, errRh, corRh, avgEta, errEta, corEta,
                    steps, params.zenoWalks, params.zenoInterval));
        } catch (Exception e) { System.err.println("Could not write results_chain.csv: " + e.getMessage()); }

        // charge.csv: P(r), same format as the star code
        double[] bins = meterZeno.getChargeBins();
        long nConf = meterZeno.getChargeConfs();
        double dr = params.chargeBinWidth;
        String chargeFile = System.getProperty("user.dir") + "/charge.csv";
        boolean chargeExists = new java.io.File(chargeFile).exists();
        try (FileWriter fw = new FileWriter(chargeFile, true)) {
            if (!chargeExists) fw.write("nSpheres,temperature,ideal,Rh,binWidth,nConf,r,r_over_Rh,C3P\n");
            for (int k = 0; k < bins.length; k++) {
                double r = (k + 0.5) * dr;
                double rho = (bins[k] / nConf) / (4 * Math.PI * r * r * dr);
                fw.write(String.format("%d,%.4f,%b,%.8f,%.4f,%d,%.6f,%.8f,%.8f%n",
                        nSpheres, temperature, ideal, avgRh, dr, nConf, r, r / avgRh, avgRh * avgRh * rho));
            }
        } catch (Exception e) { System.err.println("Could not write charge.csv: " + e.getMessage()); }

        long t2 = System.nanoTime();


        System.out.println("time: "+(t2-t1)/1e9);
    }

    public enum TruncationChoice {
        SIMPLE, SHIFT, FORCESHIFT
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public int nSpheres = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double eFENE = 0;
        public double rc = 0;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SHIFT;
        public double kBend = 0;
        public boolean ideal = false;
        public long zenoWalks = 2000;          // walks per ZENO call
        public int zenoInterval = 1000;        // MC steps between ZENO calls
        public double chargeBinWidth = 0.1;    // bin width for P(r); 0 = off

    }
}
