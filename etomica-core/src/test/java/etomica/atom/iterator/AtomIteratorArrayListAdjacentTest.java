/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.IteratorDirective;
import etomica.simulation.Simulation;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.LinkedList;

import static etomica.atom.iterator.IteratorTestAbstract.*;


/**
 * Unit test for AtomIteratorSequenceAdjacent
 * 
 * @author David Kofke
 *
 */
class AtomIteratorArrayListAdjacentTest {

    public AtomIteratorArrayListAdjacentTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }

    @BeforeEach
    public void setUp() {
        iteratorUp = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.UP);
        iteratorDn = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.DOWN);
        iteratorBoth = new AtomIteratorArrayListAdjacent(null);
    }
    
    @Test
    public void testIteration() {
        int nAtoms = 11;
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {1},11,new int[] {1});
        IAtomList atomList = sim.getBox(0).getMoleculeList().get(0).getChildList();

        //atom in middle of list
        IAtom atom = atomList.get(5);
        IAtom dnAtom = atomList.get(4);
        IAtom upAtom = atomList.get(6);
        LinkedList list = null;
        
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 2);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());
        Assertions.assertEquals(list.get(1), new AtomSetSinglet(dnAtom).toString());

        //atom at end of list        
        atom = atomList.get(nAtoms-1);
        dnAtom = atomList.get(nAtoms-2);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        //atom at beginning of list        
        atom = atomList.get(0);
        upAtom = atomList.get(1);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        Assertions.assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        //short list
        atomList = sim.getBox(0).getMoleculeList(sim.getSpecies(1)).get(0).getChildList();
        Assertions.assertEquals(atomList.size(), 1);
        atom = atomList.get(0);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);
        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);
        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);
        
    }
    
    AtomIteratorArrayListAdjacent iteratorUp, iteratorDn, iteratorBoth, iterator;
}
