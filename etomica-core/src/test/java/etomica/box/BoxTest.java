package etomica.box;

import etomica.atom.AtomType;
import etomica.box.storage.DoubleStorage;
import etomica.box.storage.Token;
import etomica.config.ConformationLinear;
import etomica.molecule.IMolecule;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class BoxTest {
    Box box;
    ISpecies s1;
    ISpecies s2;
    AtomType a1;
    AtomType a2;
    AtomType a3;

    @BeforeEach
    void setUp() {
        Simulation sim = new Simulation(Space3D.getInstance());
        a1 = AtomType.simple("A");
        a2 = AtomType.simple("B");
        a3 = AtomType.simple("C");
        s1 = new SpeciesBuilder(sim.getSpace())
                .addCount(a1, 2)
                .addCount(a2, 1)
                .withConformation(new ConformationLinear(sim.getSpace(), 1))
                .build();
        s2 = new SpeciesBuilder(sim.getSpace())
                .addCount(a3, 1)
                .withConformation(new ConformationLinear(sim.getSpace(), 1))
                .build();
        sim.addSpecies(s1);
        sim.addSpecies(s2);
        box = sim.makeBox();
    }

    @Test
    void addNewMolecule() {
        IMolecule mol = box.addNewMolecule(s1);
        assertEquals(1, box.getMoleculeList().size());
        assertEquals(3, box.getLeafList().size());
    }

    @Test
    void setNMolecules() {
        box.setNMolecules(s1, 100);
        box.setNMolecules(s2, 100);
        assertEquals(200, box.getMoleculeList().size());
        assertEquals(200, box.getMoleculeCount());
        assertEquals(400, box.getLeafList().size());
        assertEquals(400, box.getAtomCount());
        assertEquals(100, box.getNMolecules(s1));
        assertEquals(100, box.getNMolecules(s2));

        box.setNMolecules(s1, 50);
        box.setNMolecules(s2, 200);
        assertEquals(250, box.getMoleculeList().size());
        assertEquals(250, box.getMoleculeCount());
        assertEquals(350, box.getLeafList().size());
        assertEquals(350, box.getAtomCount());
    }

    @Test
    void removeMolecule() {
        box.setNMolecules(s1, 10);
        box.setNMolecules(s2, 10);

        IMolecule m1 = box.getMoleculeList().get(3);
        IMolecule m2 = box.getMoleculeList().get(14);

        box.removeMolecule(m1);
        box.removeMolecule(m2);

        assertEquals(18, box.getMoleculeList().size());
        assertEquals(36, box.getLeafList().size());

        box.addNewMolecule(s1);
        box.setNMolecules(s2, 10);

        assertEquals(20, box.getMoleculeCount());
        assertEquals(20, box.getMoleculeList().size());

        IMolecule m3 = box.getMoleculeList().get(19);
        box.removeMolecule(m3);
        assertEquals(19, box.getMoleculeCount());

        box.setNMolecules(s1, 0);
        box.setNMolecules(s2, 0);
        box.addNewMolecule(s2);
        assertEquals(1, box.getMoleculeCount());
    }

    @Test
    void storageTokens() {
        Token<DoubleStorage> tok = new Token<DoubleStorage>() {
            @Override
            public DoubleStorage createStorage(Space space) {
                return new DoubleStorage();
            }

            @Override
            public void init(int idx, DoubleStorage storage, Box box) {
                storage.create(idx).set(5.0);
            }
        };
        // Storage created before molecules added
        DoubleStorage storageA = box.getAtomStorage(tok);
        DoubleStorage storageM = box.getMolStorage(tok);

        box.setNMolecules(s1, 10);
        box.setNMolecules(s2, 10);

        for (int i = 0; i < box.getMoleculeCount(); i++) {
            assertEquals(5.0, storageM.get(i));
        }
        for (int i = 0; i < box.getAtomCount(); i++) {
            assertEquals(5.0, storageA.get(i));
        }

        Token<DoubleStorage> tok2 = new Token<DoubleStorage>() {
            @Override
            public DoubleStorage createStorage(Space space) {
                return new DoubleStorage();
            }

            @Override
            public void init(int idx, DoubleStorage storage, Box box) {
                storage.create(idx).set(11.0);
            }
        };

        DoubleStorage storageA2 = box.getAtomStorage(tok2);
        DoubleStorage storageM2 = box.getMolStorage(tok2);

        for (int i = 0; i < box.getMoleculeCount(); i++) {
            assertEquals(11.0, storageM2.get(i));
        }
        for (int i = 0; i < box.getAtomCount(); i++) {
            assertEquals(11.0, storageA2.get(i));
        }

    }
}