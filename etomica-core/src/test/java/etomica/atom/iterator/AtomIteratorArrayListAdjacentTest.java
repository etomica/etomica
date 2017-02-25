/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.LinkedList;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.ISimulation;
import etomica.atom.AtomSetSinglet;
import etomica.UnitTestUtil;


/**
 * Unit test for AtomIteratorSequenceAdjacent
 * 
 * @author David Kofke
 *
 */
public class AtomIteratorArrayListAdjacentTest extends IteratorTestAbstract {

    public AtomIteratorArrayListAdjacentTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }

    public void setUp() {
        iteratorUp = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.UP);
        iteratorDn = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.DOWN);
        iteratorBoth = new AtomIteratorArrayListAdjacent(null);
    }
    
    public void testIteration() {
        int nAtoms = 11;
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {1},11,new int[] {1});
        IAtomList atomList = sim.getBox(0).getMoleculeList().getMolecule(0).getChildList();

        //atom in middle of list
        IAtom atom = atomList.getAtom(5);
        IAtom dnAtom = atomList.getAtom(4);
        IAtom upAtom = atomList.getAtom(6);
        LinkedList list = null;
        
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 2);
        assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());
        assertEquals(list.get(1), new AtomSetSinglet(dnAtom).toString());

        //atom at end of list        
        atom = atomList.getAtom(nAtoms-1);
        dnAtom = atomList.getAtom(nAtoms-2);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(dnAtom).toString());

        //atom at beginning of list        
        atom = atomList.getAtom(0);
        upAtom = atomList.getAtom(1);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), new AtomSetSinglet(upAtom).toString());

        //short list
        atomList = sim.getBox(0).getMoleculeList(sim.getSpecies(1)).getMolecule(0).getChildList();
        assertEquals(atomList.getAtomCount(), 1);
        atom = atomList.getAtom(0);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }
    
    AtomIteratorArrayListAdjacent iteratorUp, iteratorDn, iteratorBoth, iterator;
}
