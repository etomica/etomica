package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.api.IAtomType;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.Molecule;
import etomica.atom.MoleculeSetSinglet;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.MoleculeIteratorSinglet;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */
public class MoleculeIteratorSingletTest extends MoleculeIteratorTestAbstract {

    public MoleculeIteratorSingletTest() {
        super();
    }
    
    public void setUp() {
        ISpecies species = new SpeciesSpheresHetero(Space3D.getInstance(), false, new IAtomType[0]);
        singletIterator = new MoleculeIteratorSinglet();
        testAtom1 = new Molecule(species);
        testAtom2 = new Molecule(species);
        list1 = makeTestList(new IMoleculeList[] {new MoleculeSetSinglet(testAtom1)});
        list2 = makeTestList(new IMoleculeList[] {new MoleculeSetSinglet(testAtom2)});
    }
    
    public void testIterator() {
        print("starting");
        LinkedList list = generalIteratorMethodTests(singletIterator);
        singletIterator.setMolecule(testAtom1);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list,list1);
        singletIterator.setMolecule(null);
        assertNull(singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        assertNull(singletIterator.getMolecule());
        assertTrue(list.size() == 0);
        singletIterator.setMolecule(testAtom2);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list2);
        singletIterator.setMolecule(testAtom1);
        assertEquals(testAtom1, singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list1);
        assertEquals(testAtom1, singletIterator.getMolecule());
    }
    
    private MoleculeIteratorSinglet singletIterator;
    private IMolecule testAtom1, testAtom2;
    private LinkedList list1, list2;

}
