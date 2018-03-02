package etomica;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

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
    int off = coords1dColMajor.length / 3;

    @Setup(Level.Trial)
    public void setUp() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        box = sim.makeBox();
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim, sim.getSpace());
        sim.addSpecies(species);
        box.setNMolecules(species, 1000);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);


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
        }

    }

//    @Benchmark
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

//    @Benchmark
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

    @Benchmark
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
}
