package etomica.nbr.cell;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.RandomPositionSourceDeformable;
import etomica.box.RandomPositionSourceRectangular;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * This tests NeighborCellManager functionality to iterator (via its NeighborIterator)
 * over pairs within its potential range.  Both rectangular and slanty boxes are tested
 * with box sizes (and shapes), random potential ranges, atom pair displacements.
 * Short (not requiring a lattice sum) and long (requiring a lattice sum) potential
 * ranges are tested separately.  Only the "all" iterator is tested.
 */
class NeighborCellManagerFastererTest {

    protected Simulation sim;
    protected Box boxRectangular, boxSlanty;

    @BeforeEach
    void setUp() {
        sim = new Simulation(Space3D.getInstance());
        sim.setRandom(new RandomMersenneTwister(2));
        ISpecies species = SpeciesGeneral.monatomic(sim.getSpace(), AtomType.simple("A"));
        sim.addSpecies(species);
        boxRectangular = sim.makeBox();
        boxSlanty = sim.makeBox(new BoundaryDeformablePeriodic(sim.getSpace()));
        boxRectangular.setNMolecules(species, 2);
        boxSlanty.setNMolecules(species, 2);
    }

    @Test
    void testRectangular() {
        RandomPositionSourceRectangular ps = new RandomPositionSourceRectangular(sim.getSpace(), sim.getRandom());
        ps.setBox(boxRectangular);
        IAtom atom1 = boxRectangular.getLeafList().get(0);
        IAtom atom2 = boxRectangular.getLeafList().get(1);
        Vector r1 = atom1.getPosition();
        Vector r2 = atom2.getPosition();
        Vector dr = boxRectangular.getSpace().makeVector();
        for (int n = 1; n <= 4; n++) {
            NeighborCellManagerFasterer ncmRectangular = new NeighborCellManagerFasterer(sim.getSpeciesManager(), boxRectangular, n, BondingInfo.noBonding());
            NeighborIterator iterator = ncmRectangular.makeNeighborIterator();
            long nSeenTotal = 0, nSuccessTotal = 0;
            for (int i = 0; i < 100; i++) {
                long nTotal = 0, nSeen = 0, nSuccess = 0;
                double rc = 1 + sim.getRandom().nextDouble() * 2;
                ncmRectangular.setPotentialRange(rc);
                double Ly = 20 + sim.getRandom().nextDouble() * 10;
                double Lz = 20 + sim.getRandom().nextDouble() * 10;
                r1.E(0);
                r2.E(0);
                boxRectangular.getBoundary().setBoxSize(Vector.of(20, Ly, Lz));
                ncmRectangular.init();
                for (int j = 0; j < 1000; j++) {
                    r1.E(ps.randomPosition());
                    dr.setRandomInSphere(sim.getRandom());
                    dr.TE(rc * 3);
                    boolean expected = dr.squared() < rc * rc;
                    r2.Ev1Pv2(r1, dr);
                    r2.PE(boxRectangular.getBoundary().centralImage(r2));
                    ncmRectangular.assignCellAll();
                    boolean[] seen = new boolean[]{false};
                    iterator.iterAllNeighbors(0, new NeighborIterator.NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            seen[0] = true;
                            assertEquals(rij.Mv1Squared(dr), 0, 1e-10);
                        }
                    });
                    if (expected) {
                        assertTrue(seen[0]);
                        nSuccess++;
                    } else {
                        nTotal++;
                        if (seen[0]) {
                            nSeen++;
                        }
                    }
                }
                if (false) {
                    System.out.println(nSuccess + " neighbors seen");
                    System.out.println(nSeen + "/" + nTotal + " non-neighbors seen");
                    System.out.println("ratio: " + nSeen / (double) nSuccess);
                }
                nSuccessTotal += nSuccess;
                nSeenTotal += nSeen;
            }
            if (false) {
                System.out.println(n + " final ratio: " + nSeenTotal / (double) nSuccessTotal);
            }
        }
    }

    @Test
    void testRectangularLS() {
        RandomPositionSourceRectangular ps = new RandomPositionSourceRectangular(sim.getSpace(), sim.getRandom());
        ps.setBox(boxRectangular);
        IAtom atom1 = boxRectangular.getLeafList().get(0);
        IAtom atom2 = boxRectangular.getLeafList().get(1);
        Vector r1 = atom1.getPosition();
        Vector r2 = atom2.getPosition();
        Vector dr = boxRectangular.getSpace().makeVector();
        for (int n = 1; n <= 4; n++) {
            NeighborCellManagerFasterer ncmRectangular = new NeighborCellManagerFasterer(sim.getSpeciesManager(), boxRectangular, n, BondingInfo.noBonding());
            NeighborIterator iterator = ncmRectangular.makeNeighborIterator();
            for (int i = 0; i < 100; i++) {
                double rc = (9 + sim.getRandom().nextDouble() * 10) * (1 + (n - 1) * 0.5);
                ncmRectangular.setPotentialRange(rc);
                double Lx = 20;
                double Ly = 20 + sim.getRandom().nextDouble() * 10;
                double Lz = 20 + sim.getRandom().nextDouble() * 10;
                r1.E(0);
                r2.E(0);
                boxRectangular.getBoundary().setBoxSize(Vector.of(Lx, Ly, Lz));
                ncmRectangular.init();
                for (int j = 0; j < 1000; j++) {
                    r1.E(ps.randomPosition());
                    dr.setRandomInSphere(sim.getRandom());
                    dr.TE(rc);
                    r2.Ev1Pv2(r1, dr);
                    r2.PE(boxRectangular.getBoundary().centralImage(r2));
                    ncmRectangular.assignCellAll();
                    int[] expectedCount = new int[]{0, 0};
                    Vector dri = boxRectangular.getSpace().makeVector();
                    for (int ix = -n; ix <= n; ix++) {
                        dri.setX(0, ix * Lx);
                        for (int iy = -n; iy <= n; iy++) {
                            dri.setX(1, iy * Ly);
                            for (int iz = -n; iz <= n; iz++) {
                                dri.setX(2, iz * Lz);
                                if (dri.squared() > 0 && dri.squared() < rc * rc) expectedCount[0]++; // self
                                if (dr.Mv1Squared(dri) < rc * rc) expectedCount[1]++; // i-j
                            }
                        }
                    }
                    int[] seen = new int[]{0, 0};
                    List<Vector>[] allSeen = new ArrayList[]{new ArrayList<>(), new ArrayList<>()};
                    iterator.iterAllNeighbors(0, new NeighborIterator.NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            if (rij.squared() < rc * rc) {
                                int j = jAtom.getLeafIndex();
                                seen[j]++;
                                Vector t = boxRectangular.getSpace().makeVector();
                                t.E(rij);
                                allSeen[j].add(t);
                            }
                        }
                    });
                    if (expectedCount[0] != seen[0] || expectedCount[1] != seen[1]) {
                        System.out.println("oops for " + n + " " + i + " " + j);
                        System.out.println("dr: " + dr);
                        System.out.println("self ");
                        for (Vector s : allSeen[0]) {
                            System.out.println(s);
                        }
                        System.out.println("i-j ");
                        for (Vector s : allSeen[1]) {
                            System.out.println(s);
                        }
                    }
                    assertEquals(expectedCount[0], seen[0]);
                    assertEquals(expectedCount[1], seen[1]);
                }
            }
        }
    }

    @Test
    void testDeformable() {
        RandomPositionSourceDeformable ps = new RandomPositionSourceDeformable(sim.getSpace(), sim.getRandom());
        ps.setBox(boxSlanty);
        IAtom atom1 = boxSlanty.getLeafList().get(0);
        IAtom atom2 = boxSlanty.getLeafList().get(1);
        Vector r1 = atom1.getPosition();
        Vector r2 = atom2.getPosition();
        Vector dr = boxSlanty.getSpace().makeVector();
        BoundaryDeformablePeriodic boundary = (BoundaryDeformablePeriodic) boxSlanty.getBoundary();
        for (int n = 1; n <= 4; n++) {
            NeighborCellManagerFasterer ncmSlanty = new NeighborCellManagerFasterer(sim.getSpeciesManager(), boxSlanty, n, BondingInfo.noBonding());
            NeighborIterator iterator = ncmSlanty.makeNeighborIterator();
            long nSeenTotal = 0, nSuccessTotal = 0;
            for (int i = 0; i < 100; i++) {
                long nTotal = 0, nSeen = 0, nSuccess = 0;
                double rc = 1 + sim.getRandom().nextDouble() * 2;
                ncmSlanty.setPotentialRange(rc);
                double Lyy = 20 + sim.getRandom().nextDouble() * 10;
                double Lyx = (sim.getRandom().nextDouble() - 0.5) * 10;
                double Lzz = 20 + sim.getRandom().nextDouble() * 10;
                double Lzy = (sim.getRandom().nextDouble() - 0.5) * 10;
                double Lzx = (sim.getRandom().nextDouble() - 0.5) * 10;
                boundary.setEdgeVector(0, Vector.of(20, 0, 0));
                boundary.setEdgeVector(1, Vector.of(Lyx, Lyy, 0));
                boundary.setEdgeVector(2, Vector.of(Lzx, Lzy, Lzz));
                r1.E(0);
                r2.E(0);
                ncmSlanty.init();
                for (int j = 0; j < 1000; j++) {
                    r1.E(ps.randomPosition());
//                    System.out.println("r1 "+r1);
                    dr.setRandomInSphere(sim.getRandom());
                    dr.TE(rc * 3);
//                    System.out.println("dr "+dr);
                    boolean expected = dr.squared() < rc * rc;
                    r2.Ev1Pv2(r1, dr);
//                    System.out.println("wrap "+boxSlanty.getBoundary().centralImage(r2));
                    r2.PE(boxSlanty.getBoundary().centralImage(r2));
//                    System.out.println("r2 "+r2);
                    ncmSlanty.assignCellAll();
                    boolean[] seen = new boolean[]{false};
                    iterator.iterAllNeighbors(0, new NeighborIterator.NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            seen[0] = true;
                            assertEquals(rij.Mv1Squared(dr), 0, 1e-10);
                        }
                    });
                    if (expected) {
                        assertTrue(seen[0]);
                        nSuccess++;
                    } else {
                        nTotal++;
                        if (seen[0]) {
                            nSeen++;
                        }
                    }
                }
                if (false) {
                    System.out.println(nSuccess + " neighbors seen");
                    System.out.println(nSeen + "/" + nTotal + " non-neighbors seen");
                    System.out.println("ratio: " + nSeen / (double) nSuccess);
                }
                nSuccessTotal += nSuccess;
                nSeenTotal += nSeen;
            }
            if (false) {
                System.out.println(n + " final ratio: " + nSeenTotal / (double) nSuccessTotal);
            }
        }
    }

    @Test
    void testSlantyLS() {
        RandomPositionSourceDeformable ps = new RandomPositionSourceDeformable(sim.getSpace(), sim.getRandom());
        ps.setBox(boxSlanty);
        IAtom atom1 = boxSlanty.getLeafList().get(0);
        IAtom atom2 = boxSlanty.getLeafList().get(1);
        Vector r1 = atom1.getPosition();
        Vector r2 = atom2.getPosition();
        Vector dr = boxSlanty.getSpace().makeVector();
        BoundaryDeformablePeriodic boundary = (BoundaryDeformablePeriodic) boxSlanty.getBoundary();
        for (int n = 1; n <= 4; n++) {
            NeighborCellManagerFasterer ncmSlanty = new NeighborCellManagerFasterer(sim.getSpeciesManager(), boxSlanty, n, BondingInfo.noBonding());
            NeighborIterator iterator = ncmSlanty.makeNeighborIterator();
            for (int i = 0; i < 100; i++) {
                double rc = (9 + sim.getRandom().nextDouble() * 10) * (1 + (n - 1) * 0.5);
                ncmSlanty.setPotentialRange(rc);
                double Lyy = 20 + sim.getRandom().nextDouble() * 10;
                double Lyx = (sim.getRandom().nextDouble() - 0.5) * 10;
                double Lzz = 20 + sim.getRandom().nextDouble() * 10;
                double Lzy = (sim.getRandom().nextDouble() - 0.5) * 10;
                double Lzx = (sim.getRandom().nextDouble() - 0.5) * 10;
                Vector Lx = Vector.of(20, 0, 0);
                Vector Ly = Vector.of(Lyx, Lyy, 0);
                Vector Lz = Vector.of(Lzx, Lzy, Lzz);
                boundary.setEdgeVector(0, Lx);
                boundary.setEdgeVector(1, Ly);
                boundary.setEdgeVector(2, Lz);
                r1.E(0);
                r2.E(0);
                ncmSlanty.init();
                for (int j = 0; j < 1000; j++) {
                    r1.E(ps.randomPosition());
                    dr.setRandomInSphere(sim.getRandom());
                    dr.TE(rc);
                    r2.Ev1Pv2(r1, dr);
                    r2.PE(boxSlanty.getBoundary().centralImage(r2));
                    ncmSlanty.assignCellAll();
                    int[] expectedCount = new int[]{0, 0};
                    Vector dri = boxSlanty.getSpace().makeVector();
                    for (int ix = -n; ix <= n; ix++) {
                        for (int iy = -n; iy <= n; iy++) {
                            for (int iz = -n; iz <= n; iz++) {
                                dri.E(0);
                                dri.PEa1Tv1(ix, Lx);
                                dri.PEa1Tv1(iy, Ly);
                                dri.PEa1Tv1(iz, Lz);
                                if (dri.squared() > 0 && dri.squared() < rc * rc) expectedCount[0]++; // self
                                if (dr.Mv1Squared(dri) < rc * rc) expectedCount[1]++; // i-j
                            }
                        }
                    }
                    int[] seen = new int[]{0, 0};
                    List<Vector>[] allSeen = new ArrayList[]{new ArrayList<>(), new ArrayList<>()};
                    iterator.iterAllNeighbors(0, new NeighborIterator.NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            if (rij.squared() < rc * rc) {
                                int j = jAtom.getLeafIndex();
                                seen[j]++;
                                Vector t = boxSlanty.getSpace().makeVector();
                                t.E(rij);
                                allSeen[j].add(t);
                            }
                        }
                    });
                    if (expectedCount[0] != seen[0] || expectedCount[1] != seen[1]) {
                        System.out.println("oops for " + n + " " + i + " " + j);
                        System.out.println("dr: " + dr);
                        System.out.println("self ");
                        for (Vector s : allSeen[0]) {
                            System.out.println(s);
                        }
                        System.out.println("i-j ");
                        for (Vector s : allSeen[1]) {
                            System.out.println(s);
                        }
                    }
                    assertEquals(expectedCount[0], seen[0]);
                    assertEquals(expectedCount[1], seen[1]);
                }
            }
        }
    }
}