package etomica;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.storage.VectorStorage;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

@State(Scope.Benchmark)
@Fork(1)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(time = 1, iterations = 5)
@Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
public class BenchBox {
    Box box;
    double[][] coords = new double[1000][3];
    double[] coords1d = new double[3 * 1000];
    double[] coords1dColMajor = new double[3 * 1000];
    IVecSys vecSys;
    int off = coords1dColMajor.length / 3;
    VectorStorage vecs;
    private AtomView[] atoms;
    private Vector[] posArr;

    private static class AtomView {
        Vector pos;
    }

    @Setup(Level.Trial)
    public void setUp() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesGeneral.monatomic(sim.getSpace(), AtomType.simpleFromSim(sim));
        sim.addSpecies(species);
        box = sim.makeBox();
        box.setNMolecules(species, 1000);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        vecs = new VectorStorage(Space3D.getInstance(), 1000);
        atoms = new AtomView[1000];
        posArr = new Vector[1000];

        for (int i = 0; i < box.getLeafList().size(); i++) {
            box.getLeafList().get(i).getPosition().assignTo(coords[i]);
            Vector pos = box.getLeafList().get(i).getPosition();
            coords1d[3 * i] = pos.getX(0);
            coords1d[3 * i + 1] = pos.getX(1);
            coords1d[3 * i + 2] = pos.getX(2);

            int len = box.getLeafList().size();
            coords1dColMajor[0 + i] = pos.getX(0);
            coords1dColMajor[len + i] = pos.getX(1);
            coords1dColMajor[2 * len + i] = pos.getX(2);

            vecs.get(i).E(box.getLeafList().get(i).getPosition());
            atoms[i] = new AtomView();
            atoms[i].pos = vecs.get(i);

            posArr[i] = box.getLeafList().get(i).getPosition();
        }
        vecSys = new VectorSystem2(coords1d);

    }

    @Benchmark
    public double benchNormalBox() {
        Vector dr = box.getSpace().makeVector();
        double sum = 0;
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            Vector v = atoms.get(i).getPosition();
            for (int j = 0; j < atoms.size(); j++) {
                dr.Ev1Mv2(v, atoms.get(j).getPosition());
                sum += dr.squared();
            }
        }
        return sum;
    }

    @Benchmark
    public double benchNormalBoxNoAtoms() {
        Vector dr = box.getSpace().makeVector();
        double sum = 0;
        for (int i = 0; i < posArr.length; i++) {
            Vector v = posArr[i];
            for (int j = 0; j < posArr.length; j++) {
                dr.Ev1Mv2(v, posArr[j]);
                sum += dr.squared();
            }
        }
        return sum;
    }

    @Benchmark
    public double benchNewBox() {
        VectorStorage tmp = new VectorStorage(Space3D.getInstance(), 1);
        Vector dr = tmp.get(0);
        double sum = 0;
        for (int i = 0; i < atoms.length; i++) {
            Vector v = atoms[i].pos;
            for (int j = 0; j < atoms.length; j++) {
                dr.Ev1Mv2(v, atoms[j].pos);
                sum += dr.squared();
            }
        }
        return sum;
    }

    @Benchmark
    public double benchNewBoxNoAtoms() {
        VectorStorage tmp = new VectorStorage(Space3D.getInstance(), 1);
        Vector dr = tmp.get(0);
        double sum = 0;
        for (int i = 0; i < vecs.size(); i++) {
            Vector v = vecs.get(i);
            for (int j = 0; j < vecs.size(); j++) {
                dr.Ev1Mv2(v, vecs.get(i));
                sum += dr.squared();
            }
        }
        return sum;
    }

//    @Benchmark
    public double benchCoords() {
        double sum = 0;
        for (int i = 0; i < coords.length; i++) {
            final double x = coords[i][0];
            final double y = coords[i][1];
            final double z = coords[i][2];
            for (int j = 0; j < coords.length; j++) {
                double dx = x - coords[j][0];
                double dy = y - coords[j][1];
                double dz = z - coords[j][2];
                sum += dx * dx + dy * dy + dz * dz;
            }
        }
        return sum;
    }

//    @Benchmark
    public double benchCoordsNoLoopInvariant() {
        double sum = 0;
        for (int i = 0; i < coords.length; i++) {
//            double x = coords[i][0];
//            double y = coords[i][1];
//            double z = coords[i][2];
            for (int j = 0; j < coords.length; j++) {
                double dx = coords[i][0] - coords[j][0];
                double dy = coords[i][1] - coords[j][1];
                double dz = coords[i][2] - coords[j][2];
                sum += dx * dx + dy * dy + dz * dz;
            }
        }
        return sum;

    }

    @Benchmark
    public double benchCoords1d() {
        double sum = 0;
        for (int i = 0; i < coords.length; i++) {
            final double x = coords1d[3 * i];
            final double y = coords1d[3 * i + 1];
            final double z = coords1d[3 * i + 2];
            for (int j = 0; j < coords.length; j++) {
                double dx = x - coords1d[3 * j];
                double dy = y - coords1d[3 * j + 1];
                double dz = z - coords1d[3 * j + 2];
                sum += dx * dx + dy * dy + dz * dz;
            }
        }
        return sum;
    }

//    @Benchmark
    public double benchCoords1dColMajor() {
        double sum = 0;
        for (int i = 0; i < coords.length; i++) {
            final double x = coords1dColMajor[i];
            final double y = coords1dColMajor[off + i];
            final double z = coords1dColMajor[2 * off + i];
            for (int j = 0; j < coords.length; j++) {
                double dx = x - coords1dColMajor[j];
                double dy = y - coords1dColMajor[off + j];
                double dz = z - coords1dColMajor[2 * off + j];
                sum += dx * dx + dy * dy + dz * dz;
            }
        }
        return sum;
    }

//    @Benchmark
    public double benchCoordsParallel() {
        return IntStream.range(0, coords.length).parallel().mapToDouble(i -> {
            double sum = 0;
            for (int j = 0; j < coords.length; j++) {
                double dx = coords[i][0] - coords[j][0];
                double dy = coords[i][1] - coords[j][1];
                double dz = coords[i][2] - coords[j][2];
                sum += dx * dx + dy * dy + dz * dz;
            }
            return sum;
        }).sum();
    }

    @Benchmark
    public double benchVecSys() {
        double sum = 0;
        for (int i = 0; i < vecSys.rows(); i++) {
            for (int j = 0; j < vecSys.rows(); j++) {
                Vector v = vecSys.diff(i, j);
                sum += v.squared();
            }
        }
        return sum;
    }

    public interface IVecSys {
        Vector diff(int v1, int v2);

        int rows();
    }

    public static final class VectorSystem2 implements IVecSys {
        private final double[] coords1d;
        private final int rows;

        public VectorSystem2(double[] coords) {
            coords1d = coords.clone();
            rows = coords1d.length / 3;
        }

        public Vector diff(int v1, int v2) {
            double dx = coords1d[3 * v1] - coords1d[3 * v2];
            double dy = coords1d[3 * v1 + 1] - coords1d[3 * v2 + 1];
            double dz = coords1d[3 * v1 + 2] - coords1d[3 * v2 + 2];
            return new Vector3D(dx, dy, dz);
        }

        public int rows() {
            return rows;
        }

    }
}
