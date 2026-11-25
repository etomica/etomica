/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.PotentialCompute;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.starpolymer.ConformationStarPolymerGraft;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.ClusterWeightAbs;
import etomica.virial.mcmove.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

/**
 * Samples only the target system of the star-polymer B2 calculation in {@link VirialjackCell}
 * (same model, same target-box moves), without a reference system or virial data.  It is meant for
 * timing and checking changes to the energy and neighbor code.
 * <p>
 * Energies come from cell lists (PotentialMoleculePairCell with NeighborCellManagerStar, as in
 * VirialjackCell; the default) or, with -incremental, from incremental updates
 * ({@link PotentialMoleculePairStar}).  In incremental mode ClusterSum's caching is turned off and
 * two calculations are stored (see PotentialMoleculePairStar); -clusterCaching and -numStates
 * change that, and -debugIncremental compares incremental calculations with full ones.  The box
 * is non-periodic by default; its size matters only for the cell lists.
 * <p>
 * It reports the time per MC step, running averages of the intramolecular energy, the sign of
 * f12 (the quantity B2 depends on), f12 and the core-core distance, measures of the neighbor-search
 * cost (fraction of atoms outside the box, candidate pairs visited per energy recalculation, and how
 * many of those are monomer pairs beyond rc), and the move acceptance ratios.
 * Every checkInterval steps (outside the timed region) it compares the energies with all-pairs
 * sums.  The configuration and MC step sizes can be saved (configOut) and reloaded
 * (configIn).  MC step sizes are adjusted only when starting fresh (no configIn); a run started
 * from a saved configuration keeps the saved step sizes fixed, so that with the same seed, runs
 * with different code follow the same trajectory as long as the energies agree exactly.
 */
public class VirialjackCellTarget {

    public static void main(String[] args) {
        Params params = new Params();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        final int numArms = params.numArms;
        final int armLength = params.armLength;
        final double temperature = params.temperature;
        if (armLength < 3) {
            throw new RuntimeException("arm length must be at least 3");
        }
        Space space = Space3D.getInstance();

        // --- model: identical to VirialjackCell ---
        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(space, numArms, armLength);
        conf.setSigma0(params.coreSigma);
        AtomType typeCore = AtomType.simple("Core");
        AtomType typeMono = AtomType.simple("Mono", 0);
        ISpecies species = new SpeciesBuilder(space)
                .addCount(typeCore, 1)
                .addCount(typeMono, numArms * armLength)
                .withConformation(conf)
                .build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        double coreRadius = params.coreSigma / 2.0;
        IPotential2 p2CoreCore = new P2ParticleParticle(coreRadius, params.sigmaAtom, params.epsilonAtom, params.density);
        IPotential2 p2MonoMono = P2LennardJones.makeTruncated(1, 1, new TruncationFactorySimple(params.rc));
        IPotential2 p2CoreMono = new P2ParticleSegment(coreRadius, params.sigmaAtom, params.epsilonAtom, params.density);

        // energies: cell lists (PotentialMoleculePairCell, as in VirialjackCell) or incremental
        // updates (PotentialMoleculePairStar, with -incremental)
        PotentialMoleculePairCell pCell = null;
        PotentialMoleculePairStar pStar = null;
        IPotentialMolecular pTarget;
        MoleculePairEnergy pairEnergy;
        IPotential2[][] atomPotentials;
        if (params.incremental) {
            pStar = new PotentialMoleculePairStar(sm, typeCore);
            pStar.setAtomPotential(typeCore, typeCore, p2CoreCore);
            pStar.setAtomPotential(typeMono, typeMono, p2MonoMono);
            pStar.setAtomPotential(typeCore, typeMono, p2CoreMono);
            pStar.fullInterval = params.fullInterval;
            pStar.setNumStates(params.numStates);
            pStar.debugCount = params.debugIncremental;
            pTarget = pStar;
            pairEnergy = pStar::energy;
            atomPotentials = pStar.getAtomPotentials();
        } else {
            pCell = new PotentialMoleculePairCell(sm);
            pCell.setAtomPotential(typeCore, typeCore, p2CoreCore);
            pCell.setAtomPotential(typeMono, typeMono, p2MonoMono);
            pCell.setAtomPotential(typeCore, typeMono, p2CoreMono);
            pTarget = pCell;
            pairEnergy = pCell::energy;
            atomPotentials = pCell.getAtomPotentials();
        }
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        ClusterSum targetCluster = VirialjackCell.makeB2Cluster(fTarget);
        targetCluster.setTemperature(temperature);
        // With incremental energies, ClusterSum's own caching is turned off (unless -clusterCaching)
        // so that the potential is asked for the energy after every rejected trial and always knows
        // the accepted configuration; see PotentialMoleculePairStar.
        boolean clusterCaching = !params.incremental || params.clusterCaching;
        targetCluster.setCaching(clusterCaching);

        // all intramolecular pairs are computed (bonded pairs included), as in VirialjackCell
        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };
        IPotential2 pBonding = new IPotential2() {
            public double getRange() {
                return 2;
            }

            public void u012add(double r2, double[] u012) {
            }
        };
        List<int[]> pairs = new ArrayList<>();
        for (int i = 0; i < numArms; i++) {
            pairs.add(new int[]{0, 1 + armLength * i});
            for (int j = 0; j < armLength - 1; j++) {
                pairs.add(new int[]{1 + armLength * i + j, 1 + armLength * i + j + 1});
            }
        }
        bondingInfo.setBondingPotentialPair(species, pBonding, pairs);

        // --- target box only ---
        Simulation sim = new Simulation(space, sm);
        sim.setRandom(new RandomMersenneTwister(params.seed));
        // A non-periodic boundary (the default, and what VirialjackCell gets from
        // SimulationVirialOverlap2) gives isolated molecules.  With a periodic one,
        // NeighborIteratorCellFaster uses nearest-image distances for pairs
        // involving a core, and MCMoveClusterRotateMoleculeMulti wraps atom positions, which breaks
        // molecules once atoms are more than boxLength/2 from the rotation centre (CORE_QUIRKS.md #8).
        BoxCluster box = params.periodicBox
                ? new BoxCluster(ClusterWeightAbs.makeWeightCluster(targetCluster), space, params.boxLength)
                : new BoxCluster(ClusterWeightAbs.makeWeightCluster(targetCluster), space);
        if (!params.periodicBox) {
            box.getBoundary().setBoxSize(Vector.of(params.boxLength, params.boxLength, params.boxLength));
        }
        sim.addBox(box);
        box.setNMolecules(species, 2);
        new ConfigurationCluster(space).initializeCoordinates(box);
        if (params.configIn != null) {
            new ConfigurationFile(params.configIn).initializeCoordinates(box);
            box.trialNotify();
            box.acceptNotify();
        }
        PotentialCompute pc;
        NeighborCellManagerStar cellManager = null;   // stays null for the incremental energies
        if (pStar != null) {
            pStar.setBox(box);
            pc = pStar.makePotentialComputeIntra();
        } else {
            pCell.setBox(box);
            pc = pCell.makePotentialComputeIntra();
            cellManager = new NeighborCellManagerStar(sm, box, 2, bondingInfo, typeCore);
            cellManager.setPotentialRange(params.rc);
            cellManager.init();
            pCell.setNeighborManager(cellManager);
        }
        IntegratorMC integrator = new IntegratorMC(pc, sim.getRandom(), temperature, box);
        if (pCell != null) pCell.integrator = integrator;

        // --- target moves: identical to VirialjackCell's box 1 ---
        IntArrayList[] bonding = new IntArrayList[1 + numArms * armLength];
        bonding[0] = new IntArrayList(numArms);
        for (int i = 0; i < numArms; i++) {
            bonding[0].add(1 + i * armLength);
            bonding[1 + i * armLength] = new IntArrayList(new int[]{0, 1 + i * armLength + 1});
            for (int j = 1; j < armLength - 1; j++) {
                bonding[1 + armLength * i + j] = new IntArrayList(new int[]{1 + armLength * i + j - 1, 1 + armLength * i + j + 1});
            }
            bonding[1 + (i + 1) * armLength - 1] = new IntArrayList(new int[]{1 + (i + 1) * armLength - 2});
        }
        int[] constraintMap = new int[]{0, 1, 0};
        MCMoveManager moveManager = integrator.getMoveManager();

        MCMoveClusterRotateMoleculeMulti mcMoveRotate = new MCMoveClusterRotateMoleculeMulti(sim.getRandom(), box);
        mcMoveRotate.setStepSize(Math.PI);
        mcMoveRotate.setCellManager(cellManager);
        moveManager.addMCMove(mcMoveRotate);
        MCMoveClusterMoleculeMulti mcMoveTranslate = new MCMoveClusterMoleculeMulti(sim.getRandom(), box);
        mcMoveTranslate.setCellManager(cellManager);
        moveManager.addMCMove(mcMoveTranslate);

        MCMoveClusterAngle angleInner = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, 2, (armLength + 1) / 2));
        MCMoveClusterAngle angleOuter = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, (armLength + 1) / 2 + 1, armLength));
        for (MCMoveClusterAngle m : new MCMoveClusterAngle[]{angleInner, angleOuter}) {
            m.setBox(box);
            m.setConstraintMap(constraintMap);
            m.setFixedCOM(false);
            m.setCellManager(cellManager);
            moveManager.addMCMove(m, 0.1);
        }
        ((MCMoveStepTracker) angleInner.getTracker()).setAcceptanceTarget(0.5);
        MCMoveClusterAngleMulti angleMultiInner = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, 2, (armLength + 1) / 2), 5);
        MCMoveClusterAngleMulti angleMultiOuter = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, armLength / 2 + 1, armLength), 10);
        for (MCMoveClusterAngleMulti m : new MCMoveClusterAngleMulti[]{angleMultiInner, angleMultiOuter}) {
            m.setBox(box);
            m.setConstraintMap(constraintMap);
            m.setFixedCOM(false);
            m.setCellManager(cellManager);
            moveManager.addMCMove(m, 0.4);
        }

        // optional crankshaft move on the inner part of the arms
        MCMoveCrankshaft crankshaft = null;
        if (params.crankshaftFreq > 0) {
            int aMax = params.crankshaftAMax < 0 ? armLength / 2 : params.crankshaftAMax;
            crankshaft = new MCMoveCrankshaft(pc, space, sim.getRandom(), 1, numArms, armLength, 0, aMax, params.crankshaftMaxSpan);
            crankshaft.setBox(box);
            crankshaft.setCellManager(cellManager);
            moveManager.addMCMove(crankshaft, params.crankshaftFreq);
        }

        for (MCMove m : moveManager.getMCMoves()) {
            ((MCMoveStepTracker) m.getTracker()).setNoisyAdjustment(false);
        }
        boolean adjustStepSizes = params.configIn == null;
        if (params.configIn != null) {
            readStepSizes(params.configIn + ".steps", moveManager);
        }

        System.out.println(numArms + " arms of length " + armLength + ", target system only, T = " + temperature
                + ", rc = " + params.rc + ", boxLength = " + params.boxLength + (params.periodicBox ? " (periodic)" : " (non-periodic)") + ", seed = " + params.seed);
        System.out.println((params.configIn == null ? "fresh start" : "configuration from " + params.configIn + ".pos")
                + "; step sizes " + (adjustStepSizes ? "adjusted" : "fixed")
                + "; energies " + (pStar != null ? "incremental (PotentialMoleculePairStar, full recalculation every "
                + params.fullInterval + ", " + params.numStates + " stored, " + params.debugIncremental + " debug comparisons)" : "cell lists (PotentialMoleculePairCell)")
                + (clusterCaching ? "" : "; cluster caching off"));

        // --- fresh start: relax the arms with intramolecular moves only, as VirialjackCell does ---
        if (params.configIn == null && params.relaxSteps > 0) {
            MCMoveClusterAngle[] am = {angleInner, angleOuter};
            MCMoveClusterAngleMulti[] amm = {angleMultiInner, angleMultiOuter};
            for (MCMoveClusterAngle m : am) m.skipW = true;
            for (MCMoveClusterAngleMulti m : amm) m.skipW = true;
            if (crankshaft != null) crankshaft.skipW = true;
            moveManager.setFrequency(mcMoveTranslate, 0);
            moveManager.setFrequency(mcMoveRotate, 0);
            moveManager.setEquilibrating(true);
            long tr = System.nanoTime();
            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, params.relaxSteps));
            System.out.println("relaxed arms for " + params.relaxSteps + " steps (" + (System.nanoTime() - tr) / 1e9 + " s)");
            for (MCMoveClusterAngle m : am) m.skipW = false;
            for (MCMoveClusterAngleMulti m : amm) m.skipW = false;
            if (crankshaft != null) crankshaft.skipW = false;
            moveManager.setFrequency(mcMoveTranslate, 1);
            moveManager.setFrequency(mcMoveRotate, 1);
        }
        if (params.equilibrateSteps > 0) {
            // adjust step sizes (e.g. for a move not in the saved step-size file), then run as usual
            moveManager.setEquilibrating(true);
            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, params.equilibrateSteps));
            System.out.println("adjusted step sizes for " + params.equilibrateSteps + " steps");
        }
        moveManager.setEquilibrating(adjustStepSizes);

        // --- metrics: uIntra, sign(f12), f12, core-core distance (from the values the sampling uses),
        // and the cost of the neighbor search in the current configuration: the fraction of atoms
        // outside the box (in the coarse outer cells), the candidate pairs the cell iterator visits
        // in one full energy recalculation, and how many of those are monomer pairs beyond rc.
        // The search metrics use a separate iterator over the current cells; this costs about one
        // extra recalculation per sampleInterval steps, so roughly 1% of the timed steps. ---
        IMoleculeList molecules = box.getMoleculeList();
        IAtomList leafAtoms = box.getLeafList();
        NeighborIterator countingIterator = cellManager == null ? null : cellManager.makeNeighborIterator();
        double rc2 = params.rc * params.rc, halfL = params.boxLength / 2;
        IDataSource metrics = new IDataSource() {
            final DataDoubleArray data = new DataDoubleArray(7);
            final DataDoubleArray.DataInfoDoubleArray info = new DataDoubleArray.DataInfoDoubleArray("metrics", Null.DIMENSION, new int[]{7});
            final DataTag tag = new DataTag();
            final long[] visits = new long[2];

            public IData getData() {
                double[] x = data.getData();
                x[0] = pairEnergy.energy(molecules.get(0), molecules.get(0)) + pairEnergy.energy(molecules.get(1), molecules.get(1));
                x[2] = targetCluster.value(box);
                x[1] = Math.signum(x[2]);
                x[3] = Math.sqrt(molecules.get(0).getChildList().get(0).getPosition().Mv1Squared(molecules.get(1).getChildList().get(0).getPosition()));
                int outside = 0;
                for (IAtom a : leafAtoms) {
                    Vector r = a.getPosition();
                    if (Math.abs(r.getX(0)) > halfL || Math.abs(r.getX(1)) > halfL || Math.abs(r.getX(2)) > halfL) outside++;
                }
                x[4] = outside / (double) leafAtoms.size();
                if (countingIterator == null) {
                    x[5] = x[6] = 0;   // not applicable; NaN would spoil all the averages
                    return data;
                }
                visits[0] = visits[1] = 0;
                for (int i = 0; i < leafAtoms.size(); i++) {
                    boolean iCore = leafAtoms.get(i).getType() == typeCore;
                    countingIterator.iterUpNeighbors(i, (aj, rij, k) -> {
                        visits[0]++;
                        if (!iCore && aj.getType() != typeCore && rij.squared() > rc2) visits[1]++;
                    });
                }
                x[5] = visits[0];
                x[6] = visits[1];
                return data;
            }

            public DataTag getTag() {
                return tag;
            }

            public IDataInfo getDataInfo() {
                return info;
            }
        };
        long numBlocks = 100;
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(Math.max(1, params.numSteps / params.sampleInterval / numBlocks));
        integrator.getEventManager().addListener(new DataPumpListener(metrics, acc, params.sampleInterval));

        // --- optional time series of conformational quantities that do not depend on the
        // molecule's orientation, for correlation analysis: for every arm of every molecule, the
        // distance from the core to monomers at several positions along the arm (1 = the grafted
        // monomer), and the span of 10-bond segments near the core, in the middle and at the end ---
        FileWriter seriesWriter = null;
        if (params.seriesFile != null) {
            int[] seriesPos = java.util.Arrays.stream(new int[]{5, 15, 30, 45, armLength}).filter(k -> k <= armLength).distinct().toArray();
            int[][] seriesSpans = {{1, Math.min(11, armLength)}, {Math.max(1, armLength / 2 - 5), Math.min(armLength, armLength / 2 + 5)}, {Math.max(1, armLength - 10), armLength}};
            try {
                seriesWriter = new FileWriter(params.seriesFile);
                StringBuilder header = new StringBuilder("# step");
                for (int m = 0; m < 2; m++) {
                    for (int a = 0; a < numArms; a++) {
                        for (int k : seriesPos) header.append(" m").append(m).append("a").append(a).append("_r").append(k);
                        for (int[] sp : seriesSpans) header.append(" m").append(m).append("a").append(a).append("_s").append(sp[0]).append("-").append(sp[1]);
                    }
                }
                seriesWriter.write(header + "\n");
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            FileWriter sw = seriesWriter;
            integrator.getEventManager().addListener(new IntegratorListenerAction(() -> {
                StringBuilder line = new StringBuilder(Long.toString(integrator.getStepCount()));
                for (int m = 0; m < 2; m++) {
                    IAtomList atoms = molecules.get(m).getChildList();
                    Vector core = atoms.get(0).getPosition();
                    for (int a = 0; a < numArms; a++) {
                        for (int k : seriesPos) {
                            line.append(String.format(" %.6g", Math.sqrt(atoms.get(1 + armLength * a + k - 1).getPosition().Mv1Squared(core))));
                        }
                        for (int[] sp : seriesSpans) {
                            Vector r1 = atoms.get(1 + armLength * a + sp[0] - 1).getPosition();
                            Vector r2 = atoms.get(1 + armLength * a + sp[1] - 1).getPosition();
                            line.append(String.format(" %.6g", Math.sqrt(r1.Mv1Squared(r2))));
                        }
                    }
                }
                try {
                    sw.write(line + "\n");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }, params.sampleInterval));
        }

        int[] numBadChecks = new int[1];

        // --- main run, in chunks; only the MC steps are timed ---
        long callCount0 = pCell != null ? pCell.callCount : 0, pairCount0 = pCell != null ? pCell.pairCount : 0;
        long chunk = Math.max(1, Math.min(params.checkInterval, params.numSteps));
        long stepsDone = 0;
        double tSteps = 0, tSinceReport = 0;
        long stepsSinceReport = 0;
        while (stepsDone < params.numSteps) {
            long n = Math.min(chunk, params.numSteps - stepsDone);
            long t0 = System.nanoTime();
            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, n));
            double dt = (System.nanoTime() - t0) / 1e9;
            tSteps += dt;
            tSinceReport += dt;
            stepsDone += n;
            stepsSinceReport += n;

            NeighborCellManagerStar cm = cellManager;
            PotentialMoleculePairStar ps = pStar;
            Runnable forceRecalculation = cm != null ? cm::assignCellAll : ps::forceFullNextCalculation;
            check(stepsDone, pairEnergy, forceRecalculation, atomPotentials, typeCore, box, numBadChecks);
            if (params.configOut != null && (stepsDone % params.saveInterval < n || stepsDone == params.numSteps)) {
                save(params.configOut, box, space, moveManager);
            }
            if (stepsDone % params.reportInterval < n || stepsDone == params.numSteps) {
                System.out.printf("step %d  %.1f us/step (last %d steps)  cumulative %.1f us/step%n",
                        stepsDone, 1e6 * tSinceReport / stepsSinceReport, stepsSinceReport, 1e6 * tSteps / stepsDone);
                tSinceReport = 0;
                stepsSinceReport = 0;
            }
        }

        if (seriesWriter != null) {
            try {
                seriesWriter.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        // --- final report ---
        System.out.println();
        System.out.printf("total time in MC steps: %.2f s for %d steps = %.2f us/step%n", tSteps, stepsDone, 1e6 * tSteps / stepsDone);
        IData avg = acc.getData(AccumulatorAverage.AVERAGE);
        IData err = acc.getData(AccumulatorAverage.ERROR);
        IData cor = acc.getData(AccumulatorAverage.BLOCK_CORRELATION);
        String[] names = {"uIntra", "sgn(f12)", "f12", "rCore", "outFrac", "visits", "beyondRc"};
        System.out.println("averages (sampled every " + params.sampleInterval + " steps, " + acc.getBlockCount() + " blocks):");
        for (int i = 0; i < (cellManager == null ? 5 : names.length); i++) {
            System.out.printf("  %-8s %.10g  err %.3g  cor %.3f%n", names[i], avg.getValue(i), err.getValue(i), cor.getValue(i));
        }
        double[] last = ((DataDoubleArray) metrics.getData()).getData();
        System.out.printf("final configuration: uIntra %.15g  u12 %.15g  f12 %.15g%n", last[0],
                pairEnergy.energy(molecules.get(0), molecules.get(1)), last[2]);
        System.out.println("move acceptance ratios and step sizes:");
        for (MCMove m : moveManager.getMCMoves()) {
            System.out.printf("  %-36s freq %.1f  acc %.4f  step %s%n", moveName(m, moveManager), moveManager.getFrequency(m),
                    m.getTracker().acceptanceRatio(), m instanceof MCMoveBoxStep ? "" + ((MCMoveBoxStep) m).getStepSize() : "-");
        }
        // counts include the energy checks, which are a small fraction for checkInterval >> 1
        if (pCell != null) {
            long calls = pCell.callCount - callCount0;
            System.out.printf("full energy recalculations: %.3f per step, %.0f nonzero pair terms each%n",
                    calls / (double) stepsDone, (pCell.pairCount - pairCount0) / (double) calls);
        } else {
            System.out.printf("energy calculations per step: %.3f incremental, %.4f full, %.3f restored after rejection%n",
                    pStar.numIncremental / (double) stepsDone, pStar.numFull / (double) stepsDone, pStar.numRestored / (double) stepsDone);
            long nm = pStar.numMolUnchanged + pStar.numMolRigid + pStar.numMolPartial + pStar.numMolFull;
            System.out.printf("intramolecular energies in incremental calculations: %.3f unchanged, %.3f rigid move, %.3f partial update, %.4f recomputed%n",
                    pStar.numMolUnchanged / (double) nm, pStar.numMolRigid / (double) nm, pStar.numMolPartial / (double) nm, pStar.numMolFull / (double) nm);
            long ni = pStar.numInterUnchanged + pStar.numInterPartial + pStar.numInterFull;
            System.out.printf("interaction energies in incremental calculations: %.3f unchanged, %.3f partial update, %.3f recomputed%n",
                    pStar.numInterUnchanged / (double) ni, pStar.numInterPartial / (double) ni, pStar.numInterFull / (double) ni);
            System.out.printf("moved atoms per partial update: %.1f intramolecular, %.1f interaction%n",
                    pStar.sumMovedIntra / (double) pStar.numMolPartial, pStar.sumMovedInter / (double) pStar.numInterPartial);
        }
        System.out.println("energy checks: " + (numBadChecks[0] == 0 ? "all passed" : numBadChecks[0] + " FAILED"));
    }

    /**
     * Compares the cell-based energies with all-pairs sums over the isolated molecules (no
     * periodic images): first the values the sampling is using (cached for the current
     * configuration), then values recomputed from scratch.  A failure also reports the all-pairs
     * sums with nearest-image distances for pairs involving a core, which is what
     * NeighborIteratorCellFaster does; agreement with those means the cells found every pair and
     * the difference comes from imaging alone.
     */
    private static void check(long step, MoleculePairEnergy pCell, Runnable forceRecalculation, IPotential2[][] pots,
                              AtomType typeCore, BoxCluster box, int[] numBad) {
        IMoleculeList molecules = box.getMoleculeList();
        IMolecule m0 = molecules.get(0), m1 = molecules.get(1);
        double[] raw = allPairs(box, pots, typeCore, false);
        for (int pass = 0; pass < 2; pass++) {
            if (pass == 1) {
                forceRecalculation.run();
                box.trialNotify();
                box.acceptNotify();
            }
            double u12 = pCell.energy(m0, m1);
            double uIntra = pCell.energy(m0, m0) + pCell.energy(m1, m1);
            if (differ(u12, raw[0]) || differ(uIntra, raw[1])) {
                numBad[0]++;
                double[] imaged = allPairs(box, pots, typeCore, true);
                System.out.printf("CHECK FAILED at step %d (%s): u12 cell %.15g all-pairs %.15g (core-imaged %.15g); uIntra cell %.15g all-pairs %.15g (core-imaged %.15g)%n",
                        step, pass == 0 ? "cached" : "recomputed", u12, raw[0], imaged[0], uIntra, raw[1], imaged[1]);
            }
        }
    }

    /**
     * Returns {u12, uIntra} summed over all atom pairs, with plain distances, or with
     * nearest-image distances for pairs involving a core if imageCore is true.
     */
    private static double[] allPairs(Box box, IPotential2[][] pots, AtomType typeCore, boolean imageCore) {
        IAtomList atoms = box.getLeafList();
        Vector dr = box.getSpace().makeVector();
        double u12 = 0, uIntra = 0;
        for (int i = 0; i < atoms.size(); i++) {
            IAtom a = atoms.get(i);
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom b = atoms.get(j);
                dr.Ev1Mv2(b.getPosition(), a.getPosition());
                if (imageCore && (a.getType() == typeCore || b.getType() == typeCore)) {
                    box.getBoundary().nearestImage(dr);
                }
                double u = pots[a.getType().getIndex()][b.getType().getIndex()].u(dr.squared());
                if (a.getParentGroup() == b.getParentGroup()) uIntra += u;
                else u12 += u;
            }
        }
        return new double[]{u12, uIntra};
    }

    /** energy of a pair of molecules, or the intramolecular energy when both are the same */
    private interface MoleculePairEnergy {
        double energy(IMolecule m1, IMolecule m2);
    }

    private static boolean differ(double a, double b) {
        return Math.abs(a - b) > 1e-9 * (1 + Math.abs(a) + Math.abs(b));
    }

    private static String moveName(MCMove m, MCMoveManager mm) {
        return mm.getMCMoves().indexOf(m) + ":" + m.getClass().getSimpleName();
    }

    /**
     * Writes all atom positions (unwrapped) to name.pos and the MC step sizes to name.steps,
     * each via a temporary file so an interrupted write never leaves a partial file.
     */
    private static void save(String name, Box box, Space space, MCMoveManager moveManager) {
        WriteConfiguration writer = new WriteConfiguration(space);
        writer.setBox(box);
        writer.setDoApplyPBC(false);
        writer.setFileName(name + ".pos.tmp");
        writer.actionPerformed();
        try (FileWriter fw = new FileWriter(name + ".steps.tmp")) {
            for (MCMove m : moveManager.getMCMoves()) {
                if (m instanceof MCMoveBoxStep) {
                    fw.write(moveName(m, moveManager) + " " + ((MCMoveBoxStep) m).getStepSize() + "\n");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        try {
            Files.move(Paths.get(name + ".pos.tmp"), Paths.get(name + ".pos"), StandardCopyOption.REPLACE_EXISTING);
            Files.move(Paths.get(name + ".steps.tmp"), Paths.get(name + ".steps"), StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void readStepSizes(String fileName, MCMoveManager moveManager) {
        File f = new File(fileName);
        if (!f.exists()) {
            System.out.println("no step-size file " + fileName + "; using default step sizes");
            return;
        }
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] s = line.trim().split(" +");
                for (MCMove m : moveManager.getMCMoves()) {
                    if (moveName(m, moveManager).equals(s[0])) {
                        ((MCMoveBoxStep) m).setStepSize(Double.parseDouble(s[1]));
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static class Params extends ParameterBase {
        public int numArms = 8;
        public int armLength = 60;
        public double temperature = 4;
        public double coreSigma = 1.0;
        public double sigmaAtom = 1.0;
        public double epsilonAtom = 1.0;
        public double density = 1;
        public double rc = 2.5;
        public double boxLength = 75;
        /** use incremental energies (PotentialMoleculePairStar) instead of cell lists */
        public boolean incremental = false;
        /** with incremental: number of stored calculations to start updates from (2 is enough with
         *  cluster caching off; use 4 with -clusterCaching) */
        public int numStates = 2;
        /** with incremental: keep ClusterSum's own caching of the cluster value (off by default) */
        public boolean clusterCaching = false;
        /** with incremental: full recalculation every this many calculations */
        public int fullInterval = 1000;
        /** with incremental: compare this many incremental calculations with full ones (throws on mismatch) */
        public long debugIncremental = 0;
        /** use a periodic boundary (wrong for these molecules; see the comment where the box is made) */
        public boolean periodicBox = false;
        public long numSteps = 1000000;
        public int seed = 1;
        /** steps of intramolecular-only relaxation on a fresh start */
        public long relaxSteps = 10000;
        /** metrics are sampled at this interval (steps) */
        public int sampleInterval = 100;
        /** energy checks run at this interval (steps), outside the timed region */
        public long checkInterval = 10000;
        public long reportInterval = 100000;
        public long saveInterval = 100000;
        /** read name.pos (and name.steps) to start from a saved configuration */
        public String configIn = null;
        /** write name.pos and name.steps every saveInterval steps and at the end */
        public String configOut = "jackTarget";
        /** frequency of the crankshaft move (0 = not used); the other moves total 3 */
        public double crankshaftFreq = 0;
        /** largest position along the arm for the crankshaft's first atom (0 = core); -1 = armLength/2 */
        public int crankshaftAMax = -1;
        /** largest number of bonds spanned by the crankshaft */
        public int crankshaftMaxSpan = 4;
        /** steps with step-size adjustment before the run (not timed or recorded) */
        public long equilibrateSteps = 0;
        /** write a time series of conformational quantities (see main) to this file, every sampleInterval steps */
        public String seriesFile = null;
    }
}
