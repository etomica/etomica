package etomica;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.CellLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.ApiAACell;
import etomica.nbr.cell.ApiAACellFast;
import etomica.nbr.cell.Cell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.stream.IntStream;

@State(Scope.Benchmark)
@Fork(1)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(time = 1, iterations = 5)
@Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
@SuppressWarnings("ALL") // duplicated code
public class BenchApiAACell {

    @Param("32000")
    int numAtoms;

    private ApiAACell pairIterator;
    private ApiAACellFast pairIteratorFast;
    private CellLattice lattice;
    private ExecutorService pool;

    @Setup
    public void setUp() {
        Simulation sim = new Simulation(Space3D.getInstance());
        Species s = new SpeciesSpheresMono(sim, sim.getSpace());
        sim.addSpecies(s);
        Box box = sim.makeBox();
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        box.getBoundary().setBoxSize(Vector.of(l, l, l));
        box.setNMolecules(s, numAtoms);
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(sim.getSpace()), sim.getSpace());
        configuration.initializeCoordinates(box);
        NeighborCellManager ncm = new NeighborCellManager(box, 2);
        pairIterator = new ApiAACell(3, 2, box);
        pairIterator.setLattice(ncm.getLattice());
        pairIteratorFast = new ApiAACellFast(box);
        pairIteratorFast.setLattice(ncm.getLattice());
        lattice = ncm.getLattice();
        pool = Executors.newFixedThreadPool(4);
    }

    @TearDown
    public void tearDown() {
        pool.shutdownNow();
    }

    @Benchmark
    public double benchPairIterator() {
        double sum = 0;
        pairIterator.reset();
        for(IAtomList pair = pairIterator.nextPair(); pair != null; pair = pairIterator.nextPair()) {
            sum += pair.get(0).getPosition().getX(0);
            sum += pair.get(1).getPosition().getX(0);
        }
        return sum;
    }

    @Benchmark
    public double benchPairIteratorFast() {
        double sum = 0;
        int[][] nbrCells = lattice.getUpNeighbors();
        Object[] sites = lattice.sites();
        for (int centralCellIdx = 0; centralCellIdx < sites.length; centralCellIdx++) {
            Cell centralCell = (Cell) sites[centralCellIdx];
            IAtomList centralCellAtoms = centralCell.occupants();

            if (!centralCellAtoms.isEmpty()) {

                // loop over pairs within the cell
                for (int i = 0; i < centralCellAtoms.size(); i++) {
                    for (int j = i + 1; j < centralCellAtoms.size(); j++) {
                        sum += centralCellAtoms.get(i).getPosition().getX(0);
                        sum += centralCellAtoms.get(j).getPosition().getX(0);
                    }
                }

                for (int nbrCellIdx : nbrCells[centralCellIdx]) {
                    Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
                    IAtomList nbrCellAtoms = nbrCell.occupants();

                    for (int i = 0; i < centralCellAtoms.size(); i++) {
                        for (int j = 0; j < nbrCellAtoms.size(); j++) {
                            sum += centralCellAtoms.get(i).getPosition().getX(0);
                            sum += nbrCellAtoms.get(j).getPosition().getX(0);
                        }
                    }
                }
            }
        }
        return sum;
    }

    @Benchmark
    public double benchParallelStreams() {
        int[][] nbrCells = lattice.getUpNeighbors();
        Object[] sites = lattice.sites();
        return IntStream.range(0, sites.length).parallel().mapToDouble(centralCellIdx -> {
            double sum = 0;
            Cell centralCell = (Cell) sites[centralCellIdx];
            IAtomList centralCellAtoms = centralCell.occupants();

            if (!centralCellAtoms.isEmpty()) {

                // loop over pairs within the cell
                for (int i = 0; i < centralCellAtoms.size(); i++) {
                    for (int j = i + 1; j < centralCellAtoms.size(); j++) {
                        sum += centralCellAtoms.get(i).getPosition().getX(0);
                        sum += centralCellAtoms.get(j).getPosition().getX(0);
                    }
                }

                for (int nbrCellIdx : nbrCells[centralCellIdx]) {
                    Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
                    IAtomList nbrCellAtoms = nbrCell.occupants();

                    for (int i = 0; i < centralCellAtoms.size(); i++) {
                        for (int j = 0; j < nbrCellAtoms.size(); j++) {
                            sum += centralCellAtoms.get(i).getPosition().getX(0);
                            sum += nbrCellAtoms.get(j).getPosition().getX(0);
                        }
                    }
                }
            }
            return sum;
        }).sum();
    }

    @Benchmark
    public double benchExecutorService() throws InterruptedException, ExecutionException {
        List<Callable<Double>> tasks = new ArrayList<>();
        int partitionSize = lattice.sites().length / 4;
        for (int i = 0; i < lattice.sites().length; i+=partitionSize) {
            tasks.add(new NeighborUpdater(i, Math.min(i + partitionSize, lattice.sites().length), lattice));
        }
        List<Future<Double>> results = pool.invokeAll(tasks);
        double sum = 0;
        for (Future<Double> result : results) {
            sum += result.get();
        }
        return sum;
    }

    @Benchmark
    public double benchCustomFJTask() {
        FJNeighborUpdater f = new FJNeighborUpdater(0, lattice.sites().length, lattice);
        return ForkJoinPool.commonPool().invoke(f);
    }

    private static class NeighborUpdater implements Callable<Double> {
        private final int from, to;
        private final CellLattice lattice;
        public NeighborUpdater(int from, int to, CellLattice lattice) {
            this.from = from;
            this.to = to;
            this.lattice = lattice;
        }

        @Override
        public Double call() throws Exception {
            double sum = 0;
            int[][] nbrCells = lattice.getUpNeighbors();
            Object[] sites = lattice.sites();
            for (int centralCellIdx = from; centralCellIdx < to; centralCellIdx++) {
                Cell centralCell = (Cell) sites[centralCellIdx];
                IAtomList centralCellAtoms = centralCell.occupants();

                if (!centralCellAtoms.isEmpty()) {

                    // loop over pairs within the cell
                    for (int i = 0; i < centralCellAtoms.size(); i++) {
                        for (int j = i + 1; j < centralCellAtoms.size(); j++) {
                            sum += centralCellAtoms.get(i).getPosition().getX(0);
                            sum += centralCellAtoms.get(j).getPosition().getX(0);
                        }
                    }

                    for (int nbrCellIdx : nbrCells[centralCellIdx]) {
                        Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
                        IAtomList nbrCellAtoms = nbrCell.occupants();

                        for (int i = 0; i < centralCellAtoms.size(); i++) {
                            for (int j = 0; j < nbrCellAtoms.size(); j++) {
                                sum += centralCellAtoms.get(i).getPosition().getX(0);
                                sum += nbrCellAtoms.get(j).getPosition().getX(0);
                            }
                        }
                    }
                }
            }
            return sum;
        }
    }

    private static class FJNeighborUpdater extends RecursiveTask<Double> {
        private static final int SEQ_CUTOFF = 1000;
        private final int from, to;
        private final CellLattice lattice;
        public FJNeighborUpdater(int from, int to, CellLattice lattice) {
            this.from = from;
            this.to = to;
            this.lattice = lattice;
        }

        @Override
        protected Double compute() {
            if(to - from < SEQ_CUTOFF) {
                return computeDirectly();
            } else {
                FJNeighborUpdater f1 = new FJNeighborUpdater(from, from + ((to - from) / 2), lattice);
                f1.fork();

                FJNeighborUpdater f2 = new FJNeighborUpdater(from + ((to - from) / 2), to, lattice);
                return f2.compute() + f1.join();
            }
        }

        private double computeDirectly() {
            double sum = 0;
            int[][] nbrCells = lattice.getUpNeighbors();
            Object[] sites = lattice.sites();
            for (int centralCellIdx = from; centralCellIdx < to; centralCellIdx++) {
                Cell centralCell = (Cell) sites[centralCellIdx];
                IAtomList centralCellAtoms = centralCell.occupants();

                if (!centralCellAtoms.isEmpty()) {

                    // loop over pairs within the cell
                    for (int i = 0; i < centralCellAtoms.size(); i++) {
                        for (int j = i + 1; j < centralCellAtoms.size(); j++) {
                            sum += centralCellAtoms.get(i).getPosition().getX(0);
                            sum += centralCellAtoms.get(j).getPosition().getX(0);
                        }
                    }

                    for (int nbrCellIdx : nbrCells[centralCellIdx]) {
                        Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
                        IAtomList nbrCellAtoms = nbrCell.occupants();

                        for (int i = 0; i < centralCellAtoms.size(); i++) {
                            for (int j = 0; j < nbrCellAtoms.size(); j++) {
                                sum += centralCellAtoms.get(i).getPosition().getX(0);
                                sum += nbrCellAtoms.get(j).getPosition().getX(0);
                            }
                        }
                    }
                }
            }
            return sum;
        }
    }

    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchApiAACell.class.getSimpleName())
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .build();

        new Runner(opts).run();
    }
}
