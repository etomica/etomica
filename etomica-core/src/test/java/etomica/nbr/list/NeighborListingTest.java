package etomica.nbr.list;

import etomica.action.BoxInflate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.*;

class NeighborListingTest {
    private NeighborListManager nlm;
    private PotentialMasterList pm;
    private Box box;

    private Map<IAtom, Set<IAtom>> pmNbrsSame;
    private Map<IAtom, Set<IAtom>> pmNbrsAB;

    private static Set<IAtom> getSameTypeNbrs(Box box, IAtom atom, double range) {
        Set<IAtom> sameTypeNbrs = box.getLeafList().stream()
                .filter(a -> a.getType() == atom.getType() && a != atom)
                .filter(a -> {
                    Vector dr = box.getSpace().makeVector();
                    dr.Ev1Mv2(a.getPosition(), atom.getPosition());
                    box.getBoundary().nearestImage(dr);
                    return dr.squared() < range * range;
                }).collect(Collectors.toSet());
        return sameTypeNbrs;
    }

    private static Set<IAtom> getOtherTypeNbrs(Box box, IAtom atom, double range) {
        Set<IAtom> otherTypeNbrs = box.getLeafList().stream()
                .filter(a -> a.getType() != atom.getType())
                .filter(a -> {
                    Vector dr = box.getSpace().makeVector();
                    dr.Ev1Mv2(a.getPosition(), atom.getPosition());
                    box.getBoundary().nearestImage(dr);
                    return dr.squared() < range * range;
                }).collect(Collectors.toSet());
        return otherTypeNbrs;
    }

    @BeforeEach
    void setup() {
        Simulation sim = new Simulation(Space3D.getInstance());
        SpeciesSpheresMono speciesA = new SpeciesSpheresMono(sim.getSpace(), new AtomType(new ElementSimple("A")));
        speciesA.setIsDynamic(true);
        SpeciesSpheresMono speciesB = new SpeciesSpheresMono(sim.getSpace(), new AtomType(new ElementSimple("B")));
        speciesB.setIsDynamic(true);
        sim.addSpecies(speciesA);
        sim.addSpecies(speciesB);
        box = sim.makeBox();
        box.setNMolecules(speciesA, 10);
        box.setNMolecules(speciesB, 10);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(sim.getSpace()), sim.getSpace());
        config.initializeCoordinates(box);
        pm = new PotentialMasterList(sim, 5, sim.getSpace());

        pmNbrsSame = new HashMap<>();
        pmNbrsAB = new HashMap<>();

        IPotentialAtomic mockPotentialAA = new MockPotential(2, 5) {
            @Override
            public double energy(IAtomList atoms) {
                assertEquals(atoms.get(0).getType(), speciesA.getLeafType());
                assertEquals(atoms.get(1).getType(), speciesA.getLeafType());
                assertTrue(pmNbrsSame.computeIfAbsent(atoms.get(0), a -> new HashSet<>()).add(atoms.get(1)));
                assertTrue(pmNbrsSame.computeIfAbsent(atoms.get(1), a -> new HashSet<>()).add(atoms.get(0)));
                return 0;
            }
        };

        IPotentialAtomic mockPotentialAB = new MockPotential(2, 5) {
            @Override
            public double energy(IAtomList atoms) {
                if (atoms.get(0).getType() == speciesA.getLeafType()) {
                    assertEquals(atoms.get(1).getType(), speciesB.getLeafType());
                } else {
                    assertEquals(atoms.get(1).getType(), speciesA.getLeafType());
                }
                assertTrue(pmNbrsAB.computeIfAbsent(atoms.get(0), a -> new HashSet<>()).add(atoms.get(1)));
                assertTrue(pmNbrsAB.computeIfAbsent(atoms.get(1), a -> new HashSet<>()).add(atoms.get(0)));
                return 0;
            }
        };

        IPotentialAtomic mockPotentialBB = new MockPotential(2, 5) {
            @Override
            public double energy(IAtomList atoms) {
                assertEquals(atoms.get(0).getType(), speciesB.getLeafType());
                assertEquals(atoms.get(1).getType(), speciesB.getLeafType());
                assertTrue(pmNbrsSame.computeIfAbsent(atoms.get(0), a -> new HashSet<>()).add(atoms.get(1)));
                assertTrue(pmNbrsSame.computeIfAbsent(atoms.get(1), a -> new HashSet<>()).add(atoms.get(0)));
                return 0;
            }
        };
        pm.addPotential(mockPotentialAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        pm.addPotential(mockPotentialAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
        pm.addPotential(mockPotentialBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});
        nlm = pm.getNeighborManager(box);
        nlm.reset();
    }

    @Test
    void testNeighborDistance() {
        for (IAtom atom : box.getLeafList()) {
            int typeIdx = atom.getType().getIndex();
            int otherTypeIdx = typeIdx == 1 ? 0 : 1;
            Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, 5);
            Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, 5);

            Set<IAtom> sameTypeNbrList = new HashSet<>();
            sameTypeNbrList.addAll(nlm.getUpList(atom)[atom.getType().getIndex()]);
            sameTypeNbrList.addAll(nlm.getDownList(atom)[atom.getType().getIndex()]);

            Set<IAtom> otherTypeNbrList = new HashSet<>();
            otherTypeNbrList.addAll(nlm.getUpList(atom)[otherTypeIdx]);
            otherTypeNbrList.addAll(nlm.getDownList(atom)[otherTypeIdx]);

            assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
            assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
        }
    }

    @Test
    void testPotentialMasterNbrs() {

        pm.calculate(box, new IteratorDirective(), new PotentialCalculation() {
            @Override
            public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
                potential.energy(atoms);
            }
        });

        for (IAtom atom : box.getLeafList()) {
            Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, 5);
            Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, 5);

            assertEquals(sameTypeNbrs, pmNbrsSame.get(atom), atom + " same-type neighbors");
            assertEquals(otherTypeNbrs, pmNbrsAB.get(atom), atom + " other-type neighbors");
        }
    }

    @Nested
    class AfterBoxInflate {

        @BeforeEach
        void inflateBox() {
            BoxInflate inflate = new BoxInflate(box, box.getSpace());
            inflate.setScale(1.2);
            inflate.actionPerformed();
            nlm.updateNbrsIfNeeded();
        }

        @Test
        void testNbrDistance() {
            for (IAtom atom : box.getLeafList()) {
                int typeIdx = atom.getType().getIndex();
                int otherTypeIdx = typeIdx == 1 ? 0 : 1;
                Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, 5);
                Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, 5);

                Set<IAtom> sameTypeNbrList = new HashSet<>();
                sameTypeNbrList.addAll(nlm.getUpList(atom)[atom.getType().getIndex()]);
                sameTypeNbrList.addAll(nlm.getDownList(atom)[atom.getType().getIndex()]);

                Set<IAtom> otherTypeNbrList = new HashSet<>();
                otherTypeNbrList.addAll(nlm.getUpList(atom)[otherTypeIdx]);
                otherTypeNbrList.addAll(nlm.getDownList(atom)[otherTypeIdx]);

                assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
                assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
            }
        }

        @Test
        void testPotentialMasterNbrs() {

            pm.calculate(box, new IteratorDirective(), new PotentialCalculation() {
                @Override
                public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
                    potential.energy(atoms);
                }
            });

            for (IAtom atom : box.getLeafList()) {
                Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, 5);
                Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, 5);

                assertEquals(sameTypeNbrs, pmNbrsSame.computeIfAbsent(atom, a -> new HashSet<>()), atom + " same-type neighbors");
                assertEquals(otherTypeNbrs, pmNbrsAB.computeIfAbsent(atom, a -> new HashSet<>()), atom + " other-type neighbors");
            }
        }

    }

    private static abstract class MockPotential implements IPotentialAtomic {
        private final int nBody;
        private final double range;

        public MockPotential(int nBody, double range) {
            this.nBody = nBody;
            this.range = range;
        }

        @Override
        public double getRange() {
            return range;
        }

        @Override
        public void setBox(Box box) {}

        @Override
        public int nBody() {
            return nBody;
        }
    }
}
