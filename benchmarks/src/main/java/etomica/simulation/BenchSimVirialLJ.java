package etomica.simulation;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterChainHS;
import etomica.virial.cluster.Standard;
import etomica.virial.mcmove.MCMoveClusterAtomHSChain;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.virial.simulations.VirialLJD;
import etomica.virial.wheatley.ClusterWheatleyHS;
import etomica.virial.wheatley.ClusterWheatleySoft;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@SuppressWarnings("Duplicates")
@State(Scope.Benchmark)
@Fork(1)
public class BenchSimVirialLJ {
    SimulationVirialOverlap2 sim;

    @Setup(Level.Iteration)
    public void setUp() {
        VirialLJD.VirialLJParam params = new VirialLJD.VirialLJParam();
        params.nPoints = 4;
        params.temperature = 1;
        params.numSteps = 10000000L;
        params.doChainRef = true;
        params.doHist = false;
        sim = setupVirial(params);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 3, iterations = 4)
    @Measurement(time = 10, timeUnit = TimeUnit.SECONDS, iterations = 3)
    public long integratorStep() {
        sim.integratorOS.doStep();
        return sim.integratorOS.getStepCount();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchSimVirialLJ.class.getSimpleName())
                .jvmArgs()
                .addProfiler(StackProfiler.class)
                .build();

        new Runner(opts).run();
    }

    public static SimulationVirialOverlap2 setupVirial(VirialLJD.VirialLJParam params) {
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        final double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        boolean doHist = params.doHist;
        boolean doChainRef = params.doChainRef;

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);
//        System.out.println("sigmaHSRef: "+sigmaHSRef);
//        System.out.println("B"+nPoints+"HS: "+HSBn);
//        System.out.println("Lennard Jones overlap sampling B"+nPoints+" at T="+temperature);
        Space space = Space3D.getInstance();

        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        P2LennardJones pTarget = new P2LennardJones();
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
//        if (doChainRef) System.out.println("HS Chain reference");
        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        final ClusterAbstract targetCluster = new ClusterWheatleySoft(nPoints, fTarget, 1e-12);
        targetCluster.setTemperature(temperature);

//        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");

        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")));
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints}, temperature,refCluster,targetCluster);
        sim.init();

        if (doChainRef) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.accumulators[0].setBlockSize(1);
        }

        sim.integratorOS.setNumSubSteps(1000);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        steps /= 1000;
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        //sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/20);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/10);

        System.out.println("equilibration finished");

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }


        final HistogramSimple targHist = new HistogramSimple(200, new DoubleRange(-1, 4));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2 = cPairs.getr2(i, j);
                        double r = Math.sqrt(r2);
                        if (r > 1) {
                            r = Math.log(r);
                        }
                        else {
                            r -= 1;
                        }
                        targHist.addValue(r);
                    }
                }

            }

            public void integratorInitialized(IntegratorEvent e) {}
        };

        if (doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        if (doChainRef) sim.accumulators[0].setBlockSize(1);
        for (int i=0; i<2; i++) {
//            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }
        return sim;
    }
}
