//package etomica;
//
//import etomica.atom.AtomArrayList;
//import etomica.atom.AtomType;
//import etomica.atom.IAtom;
//import etomica.atom.IAtomList;
//import etomica.atom.iterator.ApiInterArrayList;
//import etomica.atom.iterator.ApiIntraArrayList;
//import etomica.atom.iterator.AtomLeafsetIterator;
//import etomica.box.Box;
//import etomica.config.ConfigurationLattice;
//import etomica.lattice.CellLattice;
//import etomica.lattice.LatticeCubicFcc;
//import etomica.lattice.RectangularLattice;
//import etomica.nbr.cell.Cell;
//import etomica.nbr.cell.NeighborCellManager;
//import etomica.potential.IteratorDirective;
//import etomica.simulation.Simulation;
//import etomica.space.Vector;
//import etomica.space3d.Space3D;
//import etomica.species.SpeciesGeneral;
//import org.openjdk.jmh.annotations.*;
//import org.openjdk.jmh.runner.Runner;
//import org.openjdk.jmh.runner.RunnerException;
//import org.openjdk.jmh.runner.options.Options;
//import org.openjdk.jmh.runner.options.OptionsBuilder;
//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.concurrent.*;
//import java.util.function.BiConsumer;
//import java.util.stream.IntStream;
//
//@State(Scope.Benchmark)
//@Fork(1)
//@BenchmarkMode(Mode.Throughput)
//@OutputTimeUnit(TimeUnit.SECONDS)
//@Warmup(time = 1, iterations = 5)
//@Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
//@SuppressWarnings("ALL") // duplicated code
//public class BenchApiAACell {
//
//    @Param("32000")
//    int numAtoms;
//
//    private ApiAACell pairIterator;
//    private ApiAACellFast pairIteratorFast;
//    private CellLattice lattice;
//    private ExecutorService pool;
//
//    @Setup
//    public void setUp() {
//        Simulation sim = new Simulation(Space3D.getInstance());
//        SpeciesGeneral s = SpeciesGeneral.monatomic(sim.getSpace(), AtomType.simpleFromSim(sim));
//        sim.addSpecies(s);
//        Box box = sim.makeBox();
//        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
//        box.getBoundary().setBoxSize(Vector.of(l, l, l));
//        box.setNMolecules(s, numAtoms);
//        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(sim.getSpace()), sim.getSpace());
//        configuration.initializeCoordinates(box);
//        NeighborCellManager ncm = new NeighborCellManager(box, 2);
//        pairIterator = new ApiAACell(3, 2, box);
//        pairIterator.setLattice(ncm.getLattice());
//        pairIteratorFast = new ApiAACellFast(box);
//        pairIteratorFast.setLattice(ncm.getLattice());
//        lattice = ncm.getLattice();
//        pool = Executors.newFixedThreadPool(4);
//    }
//
//    @TearDown
//    public void tearDown() {
//        pool.shutdownNow();
//    }
//
//    @Benchmark
//    public double benchPairIterator() {
//        double sum = 0;
//        pairIterator.reset();
//        for(IAtomList pair = pairIterator.nextPair(); pair != null; pair = pairIterator.nextPair()) {
//            sum += pair.get(0).getPosition().getX(0);
//            sum += pair.get(1).getPosition().getX(0);
//        }
//        return sum;
//    }
//
//    @Benchmark
//    public double benchPairIteratorFast() {
//        double sum = 0;
//        int[][] nbrCells = lattice.getUpNeighbors();
//        Object[] sites = lattice.sites();
//        for (int centralCellIdx = 0; centralCellIdx < sites.length; centralCellIdx++) {
//            Cell centralCell = (Cell) sites[centralCellIdx];
//            IAtomList centralCellAtoms = centralCell.occupants();
//
//            if (!centralCellAtoms.isEmpty()) {
//
//                // loop over pairs within the cell
//                for (int i = 0; i < centralCellAtoms.size(); i++) {
//                    for (int j = i + 1; j < centralCellAtoms.size(); j++) {
//                        sum += centralCellAtoms.get(i).getPosition().getX(0);
//                        sum += centralCellAtoms.get(j).getPosition().getX(0);
//                    }
//                }
//
//                for (int nbrCellIdx : nbrCells[centralCellIdx]) {
//                    Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
//                    IAtomList nbrCellAtoms = nbrCell.occupants();
//
//                    for (int i = 0; i < centralCellAtoms.size(); i++) {
//                        for (int j = 0; j < nbrCellAtoms.size(); j++) {
//                            sum += centralCellAtoms.get(i).getPosition().getX(0);
//                            sum += nbrCellAtoms.get(j).getPosition().getX(0);
//                        }
//                    }
//                }
//            }
//        }
//        return sum;
//    }
//
//    @Benchmark
//    public double benchParallelStreams() {
//        int[][] nbrCells = lattice.getUpNeighbors();
//        Object[] sites = lattice.sites();
//        return IntStream.range(0, sites.length).parallel().mapToDouble(centralCellIdx -> {
//            double sum = 0;
//            Cell centralCell = (Cell) sites[centralCellIdx];
//            IAtomList centralCellAtoms = centralCell.occupants();
//
//            if (!centralCellAtoms.isEmpty()) {
//
//                // loop over pairs within the cell
//                for (int i = 0; i < centralCellAtoms.size(); i++) {
//                    for (int j = i + 1; j < centralCellAtoms.size(); j++) {
//                        sum += centralCellAtoms.get(i).getPosition().getX(0);
//                        sum += centralCellAtoms.get(j).getPosition().getX(0);
//                    }
//                }
//
//                for (int nbrCellIdx : nbrCells[centralCellIdx]) {
//                    Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
//                    IAtomList nbrCellAtoms = nbrCell.occupants();
//
//                    for (int i = 0; i < centralCellAtoms.size(); i++) {
//                        for (int j = 0; j < nbrCellAtoms.size(); j++) {
//                            sum += centralCellAtoms.get(i).getPosition().getX(0);
//                            sum += nbrCellAtoms.get(j).getPosition().getX(0);
//                        }
//                    }
//                }
//            }
//            return sum;
//        }).sum();
//    }
//
//    @Benchmark
//    public double benchExecutorService() throws InterruptedException, ExecutionException {
//        List<Callable<Double>> tasks = new ArrayList<>();
//        int partitionSize = lattice.sites().length / 4;
//        for (int i = 0; i < lattice.sites().length; i+=partitionSize) {
//            tasks.add(new NeighborUpdater(i, Math.min(i + partitionSize, lattice.sites().length), lattice));
//        }
//        List<Future<Double>> results = pool.invokeAll(tasks);
//        double sum = 0;
//        for (Future<Double> result : results) {
//            sum += result.get();
//        }
//        return sum;
//    }
//
//    @Benchmark
//    public double benchCustomFJTask() {
//        FJNeighborUpdater f = new FJNeighborUpdater(0, lattice.sites().length, lattice);
//        return ForkJoinPool.commonPool().invoke(f);
//    }
//
//    private static class NeighborUpdater implements Callable<Double> {
//        private final int from, to;
//        private final CellLattice lattice;
//        public NeighborUpdater(int from, int to, CellLattice lattice) {
//            this.from = from;
//            this.to = to;
//            this.lattice = lattice;
//        }
//
//        @Override
//        public Double call() throws Exception {
//            double sum = 0;
//            int[][] nbrCells = lattice.getUpNeighbors();
//            Object[] sites = lattice.sites();
//            for (int centralCellIdx = from; centralCellIdx < to; centralCellIdx++) {
//                Cell centralCell = (Cell) sites[centralCellIdx];
//                IAtomList centralCellAtoms = centralCell.occupants();
//
//                if (!centralCellAtoms.isEmpty()) {
//
//                    // loop over pairs within the cell
//                    for (int i = 0; i < centralCellAtoms.size(); i++) {
//                        for (int j = i + 1; j < centralCellAtoms.size(); j++) {
//                            sum += centralCellAtoms.get(i).getPosition().getX(0);
//                            sum += centralCellAtoms.get(j).getPosition().getX(0);
//                        }
//                    }
//
//                    for (int nbrCellIdx : nbrCells[centralCellIdx]) {
//                        Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
//                        IAtomList nbrCellAtoms = nbrCell.occupants();
//
//                        for (int i = 0; i < centralCellAtoms.size(); i++) {
//                            for (int j = 0; j < nbrCellAtoms.size(); j++) {
//                                sum += centralCellAtoms.get(i).getPosition().getX(0);
//                                sum += nbrCellAtoms.get(j).getPosition().getX(0);
//                            }
//                        }
//                    }
//                }
//            }
//            return sum;
//        }
//    }
//
//    private static class FJNeighborUpdater extends RecursiveTask<Double> {
//        private static final int SEQ_CUTOFF = 1000;
//        private final int from, to;
//        private final CellLattice lattice;
//        public FJNeighborUpdater(int from, int to, CellLattice lattice) {
//            this.from = from;
//            this.to = to;
//            this.lattice = lattice;
//        }
//
//        @Override
//        protected Double compute() {
//            if(to - from < SEQ_CUTOFF) {
//                return computeDirectly();
//            } else {
//                FJNeighborUpdater f1 = new FJNeighborUpdater(from, from + ((to - from) / 2), lattice);
//                f1.fork();
//
//                FJNeighborUpdater f2 = new FJNeighborUpdater(from + ((to - from) / 2), to, lattice);
//                return f2.compute() + f1.join();
//            }
//        }
//
//        private double computeDirectly() {
//            double sum = 0;
//            int[][] nbrCells = lattice.getUpNeighbors();
//            Object[] sites = lattice.sites();
//            for (int centralCellIdx = from; centralCellIdx < to; centralCellIdx++) {
//                Cell centralCell = (Cell) sites[centralCellIdx];
//                IAtomList centralCellAtoms = centralCell.occupants();
//
//                if (!centralCellAtoms.isEmpty()) {
//
//                    // loop over pairs within the cell
//                    for (int i = 0; i < centralCellAtoms.size(); i++) {
//                        for (int j = i + 1; j < centralCellAtoms.size(); j++) {
//                            sum += centralCellAtoms.get(i).getPosition().getX(0);
//                            sum += centralCellAtoms.get(j).getPosition().getX(0);
//                        }
//                    }
//
//                    for (int nbrCellIdx : nbrCells[centralCellIdx]) {
//                        Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
//                        IAtomList nbrCellAtoms = nbrCell.occupants();
//
//                        for (int i = 0; i < centralCellAtoms.size(); i++) {
//                            for (int j = 0; j < nbrCellAtoms.size(); j++) {
//                                sum += centralCellAtoms.get(i).getPosition().getX(0);
//                                sum += nbrCellAtoms.get(j).getPosition().getX(0);
//                            }
//                        }
//                    }
//                }
//            }
//            return sum;
//        }
//    }
//
//    public static void main(String[] args) throws RunnerException {
//
//        Options opts = new OptionsBuilder()
//                .include(BenchApiAACell.class.getSimpleName())
////                .jvmArgs(
////                        "-XX:+UnlockDiagnosticVMOptions",
////                        "-XX:+PrintAssembly",
////                        "-XX:PrintAssemblyOptions=intel",
////                        "-XX:CompileCommand=print,*BoxBench.bench*")
//                .build();
//
//        new Runner(opts).run();
//    }
//
//    private static class ApiAACellFast {
//        private final Box box;
//        private CellLattice lattice;
//
//        public ApiAACellFast(Box box) {
//            this.box = box;
//        }
//
//        public void setLattice(CellLattice lattice) {
//            this.lattice = lattice;
//        }
//
//        public void forEachPair(AtomPairConsumer action) {
//            int[][] nbrCells = lattice.getUpNeighbors();
//            for (int centralCellIdx = 0; centralCellIdx < lattice.sites().length; centralCellIdx++) {
//
//                Cell centralCell = (Cell) lattice.sites()[centralCellIdx];
//                IAtomList centralCellAtoms = centralCell.occupants();
//
//                if (!centralCellAtoms.isEmpty()) {
//
//                    // loop over pairs within the cell
//                    for (int i = 0; i < centralCellAtoms.size(); i++) {
//                        for (int j = i + 1; j < centralCellAtoms.size(); j++) {
//                            action.accept(centralCellAtoms.get(i), centralCellAtoms.get(j));
//                        }
//                    }
//
//                    for (int nbrCellIdx : nbrCells[centralCellIdx]) {
//                        Cell nbrCell = (Cell) lattice.sites()[nbrCellIdx];
//                        IAtomList nbrCellAtoms = nbrCell.occupants();
//
//                        for (int i = 0; i < centralCellAtoms.size(); i++) {
//                            for (int j = 0; j < nbrCellAtoms.size(); j++) {
//                                action.accept(centralCellAtoms.get(i), nbrCellAtoms.get(j));
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//
//        @FunctionalInterface
//        public interface AtomPairConsumer extends BiConsumer<IAtom, IAtom> {
//            @Override
//            void accept(IAtom atom1, IAtom atom2);
//        }
//    }
//
//    /**
//     * Returns iterates formed from all cell-based neighbor pairs.
//     */
//    private static class ApiAACell {
//
//        protected final boolean[] periodicity;
//        private final Box box;
//        private final ApiIntraArrayList intraListIterator;
//        private final ApiInterArrayList interListIterator;
//        private final CellLattice.NeighborIterator neighborIterator;
//        private final RectangularLattice.Iterator cellIterator;
//        private AtomLeafsetIterator listIterator;
//
//        /**
//         * Constructor makes iterator that must have box specified and then be
//         * reset() before iteration.
//         *
//         * @param D     the dimension of the space of the simulation (used to
//         *              construct cell iterators)
//         * @param range the distance within which pairs of atoms are considered
//         *              neighbors. Used to define neighbor cells; some iterates may
//         *              exceed this separation
//         */
//        public ApiAACell(int D, double range, Box box) {
//            cellIterator = new RectangularLattice.Iterator(D);
//            neighborIterator = new CellLattice.NeighborIterator(D, range);
//            neighborIterator.setDirection(IteratorDirective.Direction.UP);
//            interListIterator = new ApiInterArrayList(new AtomArrayList(), new AtomArrayList());
//            intraListIterator = new ApiIntraArrayList();
//            listIterator = intraListIterator;
//            periodicity = new boolean[D];
//            this.box = box;
//        }
//
//        public void setLattice(CellLattice lattice) {
//            cellIterator.setLattice(lattice);
//            neighborIterator.setLattice(lattice);
//            unset();
//        }
//
//        /**
//         * Returns the number of atom pairs the iterator will return if reset and
//         * iterated in its present state.
//         */
//        public int size() {
//            int count = 0;
//            reset();
//            for (Object a = nextPair(); a != null; a = nextPair()) {
//                count++;
//            }
//            return count;
//        }
//
//        public IAtomList nextPair() {
//            IAtomList nextPair = listIterator.next();
//            if (nextPair == null) {
//                return advanceLists();
//            }
//            return nextPair;
//        }
//
//        public void unset() {
//            listIterator.unset();
//        }
//
//        /**
//         * Returns 2, indicating that this is a pair iterator.
//         */
//        public int nBody() {
//            return 2;
//        }
//
//        public void reset() {
//            for (int i = 0; i < periodicity.length; i++) {
//                periodicity[i] = box.getBoundary().getPeriodicity(i);
//            }
//            neighborIterator.setPeriodicity(periodicity);
//            cellIterator.reset();
//            neighborIterator.checkDimensions();
//            neighborIterator.unset();
//            listIterator.unset();
//            //System.out.println("reset in ApiAACell");
//        }//end of reset
//
//        // Moves to next pair of lists that can provide an iterate
//        // This should be invoked only if listIterator.hasNext is false
//        private IAtomList advanceLists() {
//            do {
//                //advance neighbor cell
//                if (neighborIterator.hasNext()) {
//                    interListIterator.setInnerList(((Cell) neighborIterator.next()).occupants());
//                    listIterator = interListIterator;
//                    interListIterator.reset();
//                    IAtomList pair = listIterator.next();
//                    if (pair != null) {
//                        return pair;
//                    }
//
//                    //advance central cell and set up neighbor cell iterator if
//                    // central cell has some molecules
//                } else if (cellIterator.hasNext()) {
//                    AtomArrayList list = ((Cell) cellIterator.peek()).occupants();
//                    neighborIterator.setSite(cellIterator.nextIndex());
//
//                    if (!list.isEmpty()) {//central cell has molecules
//                        interListIterator.setOuterList(list); //for neighbor-cell looping
//                        intraListIterator.setList(list);//for intra-cell looping
//                        neighborIterator.reset();
//
//                        listIterator = intraListIterator;
//                        intraListIterator.reset();
//
//                        IAtomList pair = listIterator.next();
//                        if (pair != null) {
//                            return pair;
//                        }
//                    } else {//no molecules in central cell
//                        neighborIterator.unset();
//                        listIterator.unset();
//                    }
//                } else {//no more cells at all
//                    return null;
//                }
//            } while (true);
//        }//end of advanceCell
//
//        /**
//         * @return Returns the cellIterator.
//         */
//        public CellLattice.NeighborIterator getNbrCellIterator() {
//            return neighborIterator;
//        }
//    }
//}
