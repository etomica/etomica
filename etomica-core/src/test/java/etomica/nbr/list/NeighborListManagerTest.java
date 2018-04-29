package etomica.nbr.list;

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
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.*;

class NeighborListManagerTest {
    private NeighborListManager nlm;
    private PotentialMasterList pm;
    private Box box;

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
        IPotentialAtomic mockPotentialAA = new IPotentialAtomic() {
            @Override
            public double energy(IAtomList atoms) {
                assertEquals(atoms.get(0).getType(), speciesA.getLeafType());
                assertEquals(atoms.get(1).getType(), speciesA.getLeafType());
                return 0;
            }

            @Override
            public double getRange() {
                return 5;
            }

            @Override
            public void setBox(Box box) {

            }

            @Override
            public int nBody() {
                return 2;
            }
        };

        IPotentialAtomic mockPotentialAB = new IPotentialAtomic() {
            @Override
            public double energy(IAtomList atoms) {
                if (atoms.get(0).getType() == speciesA.getLeafType()) {
                    assertEquals(atoms.get(1).getType(), speciesB.getLeafType());
                } else {
                    assertEquals(atoms.get(1).getType(), speciesA.getLeafType());
                }
                return 0;
            }

            @Override
            public double getRange() {
                return 5;
            }

            @Override
            public void setBox(Box box) {

            }

            @Override
            public int nBody() {
                return 2;
            }
        };
        pm.addPotential(mockPotentialAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        pm.addPotential(mockPotentialAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
        nlm = pm.getNeighborManager(box);
        nlm.reset();
    }

    @Test
    void testNeighborDistance() {
        for (IAtom atom : box.getLeafList()) {
            int typeIdx = atom.getType().getIndex();
            int otherTypeIdx = typeIdx == 1 ? 0 : 1;
            Set<IAtom> sameTypeNbrs = box.getLeafList().stream()
                    .filter(a -> a.getType() == atom.getType() && a != atom)
                    .filter(a -> {
                        Vector dr = box.getSpace().makeVector();
                        dr.Ev1Mv2(a.getPosition(), atom.getPosition());
                        box.getBoundary().nearestImage(dr);
                        return dr.squared() < 25;
                    }).collect(Collectors.toSet());

            Set<IAtom> otherTypeNbrs = box.getLeafList().stream()
                    .filter(a -> a.getType() != atom.getType())
                    .filter(a -> {
                        Vector dr = box.getSpace().makeVector();
                        dr.Ev1Mv2(a.getPosition(), atom.getPosition());
                        box.getBoundary().nearestImage(dr);
                        return dr.squared() < 25;
                    }).collect(Collectors.toSet());

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
    void testNeighborTypes() {

        pm.calculate(box, new IteratorDirective(null), new PotentialCalculation() {
            @Override
            public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
                potential.energy(atoms);
            }
        });
    }

}