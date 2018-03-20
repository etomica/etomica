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

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(time = 1, iterations = 5)
@Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
public class BenchApiAACell {

    @Param("32000")
    int numAtoms;

    private ApiAACell pairIterator;
    private ApiAACellFast pairIteratorFast;
    private CellLattice lattice;

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
