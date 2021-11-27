package etomica.nbr.list;

import etomica.action.BoxInflate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.compute.NeighborIterator;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;

class NeighborListingTest {
    private NeighborListManager nlm;
    private Box box;

    private static final double NBR_RANGE = 3.5;
    private static final double POTENTIAL_RANGE = 3.;
    private SpeciesGeneral speciesA;
    private SpeciesGeneral speciesB;
    private Simulation sim;

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
        sim = new Simulation(Space3D.getInstance());
        speciesA = SpeciesGeneral.monatomic(sim.getSpace(), new AtomType(new ElementSimple("A")), true);
        speciesB = SpeciesGeneral.monatomic(sim.getSpace(), new AtomType(new ElementSimple("B")), true);
        sim.addSpecies(speciesA);
        sim.addSpecies(speciesB);
        box = sim.makeBox();
        box.setNMolecules(speciesA, 20);
        box.setNMolecules(speciesB, 20);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(sim.getSpace()), sim.getSpace());
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones();
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(p2lj, POTENTIAL_RANGE);
        IPotential2[][] potentials = new IPotential2[][]{{p2,p2},{p2,p2}};
        nlm = new NeighborListManager(sim.getSpeciesManager(), box, 2, NBR_RANGE, BondingInfo.noBonding());
        nlm.setPairPotentials(potentials);
        nlm.setDoDownNeighbors(true);
        nlm.init();
        nlm.reset();
    }

    @Test
    void testNeighborDistance() {
        for (IAtom atom : box.getLeafList()) {
            Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, NBR_RANGE);
            Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, NBR_RANGE);

            Set<IAtom> sameTypeNbrList = new HashSet<>();
            Set<IAtom> otherTypeNbrList = new HashSet<>();
            nlm.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    if (jAtom.getType() == atom.getType()) sameTypeNbrList.add(jAtom);
                    else otherTypeNbrList.add(jAtom);
                }
            });

            assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
            assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
        }
    }

    @Nested
    class AfterBoxInflate {

        @BeforeEach
        void inflateBox() {
            BoxInflate inflate = new BoxInflate(box, box.getSpace());
            inflate.setScale(1.2);
            inflate.actionPerformed();
            nlm.init();
            nlm.reset();
        }

        @Test
        void testNbrDistance() {
            for (IAtom atom : box.getLeafList()) {
                Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, NBR_RANGE);
                Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, NBR_RANGE);

                Set<IAtom> sameTypeNbrList = new HashSet<>();
                Set<IAtom> otherTypeNbrList = new HashSet<>();
                nlm.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        if (jAtom.getType() == atom.getType()) sameTypeNbrList.add(jAtom);
                        else otherTypeNbrList.add(jAtom);
                    }
                });

                assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
                assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
            }
        }
    }
    
    @Nested
    class AfterAddAtoms {
        @BeforeEach
        void addAtoms() {
            RandomPositionSource rand = new RandomPositionSourceRectangular(box.getSpace(), sim.getRandom());
            rand.setBox(box);
            for (int i = 0; i < 50; i++) {
                box.addNewMolecule(speciesA).getChildList().get(0).getPosition().E(rand.randomPosition());
            }
            nlm.reset();
        }

        @Test
        void testNbrDistance() {
            for (IAtom atom : box.getLeafList()) {
                Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, NBR_RANGE);
                Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, NBR_RANGE);

                Set<IAtom> sameTypeNbrList = new HashSet<>();
                Set<IAtom> otherTypeNbrList = new HashSet<>();
                nlm.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        if (jAtom.getType() == atom.getType()) sameTypeNbrList.add(jAtom);
                        else otherTypeNbrList.add(jAtom);
                    }
                });

                assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
                assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
            }
        }
    }

    @Nested
    class AfterSetRange {
        @BeforeEach
        void setRange() {
            nlm.setNeighborRange(NBR_RANGE + 1);
            nlm.reset();
        }

        @Test
        void testNbrDistance() {
            for (IAtom atom : box.getLeafList()) {
                Set<IAtom> sameTypeNbrs = getSameTypeNbrs(box, atom, NBR_RANGE + 1);
                Set<IAtom> otherTypeNbrs = getOtherTypeNbrs(box, atom, NBR_RANGE + 1);

                Set<IAtom> sameTypeNbrList = new HashSet<>();
                Set<IAtom> otherTypeNbrList = new HashSet<>();
                nlm.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij) {
                        if (jAtom.getType() == atom.getType()) sameTypeNbrList.add(jAtom);
                        else otherTypeNbrList.add(jAtom);
                    }
                });

                assertEquals(sameTypeNbrs, sameTypeNbrList, atom + " same-type neighbors");
                assertEquals(otherTypeNbrs, otherTypeNbrList, atom + " other-type neighbors");
            }
        }
    }
}
