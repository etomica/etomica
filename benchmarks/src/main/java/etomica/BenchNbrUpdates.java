package etomica;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchNbrUpdates {
    public NeighborListManager nlm;
    public IAtom a;

    @Setup(Level.Iteration)
    public void setUp() {
        Simulation sim = new Simulation(Space3D.getInstance());
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim, sim.getSpace());
        sim.addSpecies(species);
        Box box = sim.makeBox(new BoundaryRectangularPeriodic(sim.getSpace(), 10000));
        box.setNMolecules(species, 100_000);
        Configuration initialConfig = new ConfigurationLattice(new LatticeCubicFcc(sim.getSpace()), sim.getSpace());
        initialConfig.initializeCoordinates(box);
        PotentialMasterList pm = new PotentialMasterList(sim, sim.getSpace());
        pm.addPotential(new P2SoftSphericalTruncated(sim.getSpace(), new P2LennardJones(sim.getSpace()), 3), new AtomType[]{species.getLeafType(), species.getLeafType()});
        pm.reset();
        RandomPositionSource rand = new RandomPositionSourceRectangular(sim.getSpace(), sim.getRandom());
        rand.setBox(box);
        Configuration randConfig = b -> {
            for (IAtom atom : b.getLeafList()) {
                atom.getPosition().E(rand.randomPosition());
            }
        };
        randConfig.initializeCoordinates(box);
        nlm = pm.getNeighborManager(box);
        a = box.getLeafList().get(0);
    }

    @Benchmark
    @BenchmarkMode(Mode.SingleShotTime)
    @OutputTimeUnit(TimeUnit.MILLISECONDS)
    @Warmup(iterations = 5)
    @Measurement(iterations = 10)
    public int benchNbrUpdates() {
        nlm.updateNbrsIfNeeded();
        return nlm.getUpList(a)[0].size();
    }

    @Benchmark
    @BenchmarkMode(Mode.SingleShotTime)
    @OutputTimeUnit(TimeUnit.MILLISECONDS)
    @Warmup(iterations = 5)
    @Measurement(iterations = 10)
    @Fork(value = 1, jvmArgsAppend = "-Detomica.nbr.parallel=true")
    public int benchNbrUpdatesParallel() {
        nlm.updateNbrsIfNeeded();
        return nlm.getUpList(a)[0].size();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchNbrUpdates.class.getSimpleName())
                .build();

        new Runner(opts).run();
    }
}
