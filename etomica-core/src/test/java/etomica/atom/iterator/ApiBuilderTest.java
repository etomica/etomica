/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.*;
import etomica.box.Box;
import etomica.simulation.Simulation;

import java.util.LinkedList;

/**
 * Tests the iterators made by the various static methods in ApiBuilder.
 * 
 * @author David Kofke
 *  
 */
public class ApiBuilderTest extends IteratorTestAbstract {

    private static final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private static final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;
    int n0a, nAtoms, n1a;
    private Simulation sim;
    private IMolecule parent;
    private IAtom target;
    private IAtom targetFirst;

    //**********  adjacent/nonadjacent setup -- basis is a leaf atom
//    private void setup2() {
//        IMoleculeList moleculeList = sim.getBox(0).getMoleculeList(sim.getSpeciesManager().getSpecies(1));
//        parent = moleculeList.getMolecule(5);//leaf-atom basis
//        target = parent;//atom5
//        targetFirst = moleculeList.getMolecule(0);//atom0 
//        targetLast = moleculeList.getMolecule(9);//atom9
//        up = moleculeList.getMolecule(6);
//        upNon = new IMolecule[] {moleculeList.getMolecule(7),moleculeList.getMolecule(8),
//                moleculeList.getMolecule(9)};
//        upFirst = moleculeList.getMolecule(1);
//        upFirstNon = new IMolecule[] {moleculeList.getMolecule(2),moleculeList.getMolecule(3),
//                moleculeList.getMolecule(4),moleculeList.getMolecule(5),moleculeList.getMolecule(6),
//                moleculeList.getMolecule(7),moleculeList.getMolecule(8),moleculeList.getMolecule(9)};
//        dn = moleculeList.getMolecule(4);
//        dnNon = new IMolecule[] {moleculeList.getMolecule(3),moleculeList.getMolecule(2),
//                moleculeList.getMolecule(1),moleculeList.getMolecule(0)};
//        dnLast = moleculeList.getMolecule(8);
//        dnLastNon = new IMolecule[] {moleculeList.getMolecule(7),moleculeList.getMolecule(6),
//                moleculeList.getMolecule(5),moleculeList.getMolecule(4),moleculeList.getMolecule(3),
//                moleculeList.getMolecule(2),moleculeList.getMolecule(1),moleculeList.getMolecule(0)};
//    }
    private IAtom targetLast;
    private IAtom up;
    private IAtom[] upNon;
    private IAtom upFirst;
    private IAtom[] upFirstNon;
    private IAtom dn;
    private IAtom[] dnNon;
    private IAtom dnLast;
    private IAtom[] dnLastNon;
    private IAtom iterate;
    private IAtom iterateFirst;
    private IAtom iterateLast;
   public ApiBuilderTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }

   /**
    * Setup used by adjacentPairIterator and nonAdjacentPairIterator tests.
    *
    */
    public void setUpA() {
        n0a = 5;
        nAtoms = 10;
        n1a = 10;
        sim = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a, 1 },
                nAtoms, new int[] { n1a, 2});
    }

    /**
     * Sets up various choices of basis, target, direction and checks that
     * non-adjacent atoms are given in pairs with targeted iterate.
     */
    public void testNonAdjacentPairIterator() {
        setUpA();
        ApiIntragroup api = ApiBuilder.makeNonAdjacentPairIterator();

        setup1();
        nonAdjacentPairTests(api);

        //XXX depends on https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=245#c0
//        setup2();
//        api.setBasis(new AtomSetSinglet(parent));
//
//        //target matches basis, no direction, two iterates
//        api.setDirection(null);
//        api.setTarget(target);
//        if (false) testApiIterates(api, target, upNon, dnNon);
//
//        //target matches basis, up, one iterate
//        api.setDirection(UP);
//        testApiIterates(api, UP, target, upNon);
//
//        //target matches basis, down, one iterate
//        api.setDirection(DOWN);
//        testApiIterates(api, DOWN, target, dnNon);
//
//        //target doesn't match basis, no iterates
//        api.setTarget(targetFirst);
//        testNoIterates(api);

        setup4();
        nonAdjacentPairTests(api);

    }

    public void testIntergroupTypeIterator() {
        //make tree of two species
        //species 0 has 5 molecules, each with 5 atoms, 3 of one type, 2 of another
        //species 1 has 7 molecules, each with 11 atoms, 4 of one type, 1 of another, and 6 of another
        //iterator must loop over pairs formed from molecules of each species
        sim = UnitTestUtil.makeMultitypeSpeciesTree(new int[]{5, 7},
                new int[][] {{3,2},{4,1,6}});
        ISpecies species0 = sim.getSpecies(0);
        ISpecies species1 = sim.getSpecies(1);
        AtomType[] types = new AtomType[2];
        MoleculePair basisPair = new MoleculePair();

        Box box = sim.getBox(0);
        IMoleculeList moleculeList0 = box.getMoleculeList(species0);
        IMoleculeList moleculeList1 = box.getMoleculeList(species1);

        //test 3-atom type and 4-atom type, no target
        basisPair.atom0 = moleculeList0.getMolecule(2);
        basisPair.atom1 = moleculeList1.getMolecule(1);
        types[0] = species0.getAtomType(0);
        types[1] = species1.getAtomType(0);
        ApiIntergroup api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 12);
        //test 3 and 4, one of the 3 given as target
        IAtom target0 = moleculeList0.getMolecule(2).getChildList().getAtom(1);
        api.setTarget(target0);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 4);
        //test 3 and 4, one of the 4 given as target
        IAtom target1 = moleculeList1.getMolecule(1).getChildList().getAtom(0);
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);
        //give target that isn't the specified type
        target0 = moleculeList0.getMolecule(2).getChildList().getAtom(4);
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = moleculeList1.getMolecule(1).getChildList().getAtom(10);
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        //same tests, but switch order of basis; nothing should give iterates
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = moleculeList0.getMolecule(2);
        basisPair.atom0 = moleculeList1.getMolecule(1);
        types[0] = species0.getAtomType(0);
        types[1] = species1.getAtomType(0);
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        testNoIterates(api);
        //test 3 and 4, one of the 3 given as target
        target0 = moleculeList0.getMolecule(2).getChildList().getAtom(1);
        api.setTarget(target0);
        testNoIterates(api);
        //test 3 and 4, one of the 4 given as target
        target1 = moleculeList1.getMolecule(1).getChildList().getAtom(0);
        api.setTarget(target1);
        testNoIterates(api);

        //same tests, but switch order of basis and switch order of types
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = moleculeList0.getMolecule(2);
        basisPair.atom0 = moleculeList1.getMolecule(1);
        types[1] = species0.getAtomType(0);
        types[0] = species1.getAtomType(0);
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(12, list0.size());
        //test 3 and 4, one of the 3 given as target
        target0 = moleculeList0.getMolecule(2).getChildList().getAtom(1);
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(4, list1.size());
        //test 3 and 4, one of the 4 given as target
        target1 = moleculeList1.getMolecule(1).getChildList().getAtom(0);
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(3, list1.size());

        //test 3-atom type and 1-atom type, no target
        basisPair.atom0 = moleculeList0.getMolecule(2);
        basisPair.atom1 = moleculeList1.getMolecule(1);
        types[0] = species0.getAtomType(0);
        types[1] = species1.getAtomType(1);
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(3, list0.size());
        //test 3 and 1, one of the 3 given as target
        target0 = moleculeList0.getMolecule(2).getChildList().getAtom(1);
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(1, list1.size());
        //test 3 and 1, the 1 given as target
        target1 = moleculeList1.getMolecule(1).getChildList().getAtom(4);
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(3, list1.size());
        //give target that isn't the specified type
        target0 = moleculeList0.getMolecule(2).getChildList().getAtom(4);
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = moleculeList1.getMolecule(1).getChildList().getAtom(10);
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        basisPair.atom0 = moleculeList0.getMolecule(2);
        basisPair.atom1 = moleculeList1.getMolecule(1);
        types[0] = species0.getAtomType(0);
        types[1] = species1.getAtomType(0);
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list1 = generalIteratorMethodTests(api);
        assertEquals(12, list1.size());

        //incomplete
    }

    public void testIntragroupTypeIterator() {
        //AtomsetIteratorBasisDependent makeIntragroupTypeIterator(AtomType[] types)
        //incomplete
    }

    public void testAdjacentPairIterator() {
        setUpA();
        ApiIntragroup api = ApiBuilder.makeAdjacentPairIterator();

        setup1();
        adjacentPairTests(api);

        //XXX depends on https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=245#c0
        //************ test basis is leaf
//        setup2();
//        api.setBasis(new AtomSetSinglet(parent));
//
//        //target matches basis, no direction, two iterates
//        api.setDirection(null);
//        api.setTarget(target);
//        testApiTwoIterates(api, new AtomPair(target, up), new AtomPair(dn, target));
//
//        //target matches basis, up, one iterate
//        api.setDirection(UP);
//        testApiOneIterate(api, new AtomPair(target, up));
//
//        //target matches basis, down, one iterate
//        api.setDirection(DOWN);
//        testApiOneIterate(api, new AtomPair(dn, target));
//
//        //target doesn't match basis, no iterates
//        api.setTarget(targetFirst);
//        testNoIterates(api);

        setup4();
        adjacentPairTests(api);

    }

    //******* adjacent/nonadjacent setup -- basis has only one child
    private void setup4() {
        Box box = sim.getBox(1);
        ISpecies species1 = sim.getSpecies(1);
        parent = box.getMoleculeList(species1).getMolecule(0);//box1, species1, molecule0
        target = parent.getChildList().getAtom(0);
        targetFirst = target;
        targetLast = target;
        up = dn = upFirst = dnLast = null;
        upNon = upFirstNon = dnNon = dnLastNon = new IAtom[0];
        iterate = target;
        iterateFirst = target;
        iterateLast = target;
    }

    //******* adjacent/nonadjacent setup -- basis has child atoms, target is among them
    private void setup1() {
        parent = sim.getBox(0).getMoleculeList(sim.getSpecies(0)).getMolecule(2);
        IAtomList childList = parent.getChildList();
        target = childList.getAtom(5);
        targetFirst = childList.getAtom(0);
        targetLast = childList.getAtom(9);
        up = childList.getAtom(6);
        upNon = new IAtom[] {childList.getAtom(7),childList.getAtom(8),childList.getAtom(9)};
        upFirst = childList.getAtom(1);
        upFirstNon = new IAtom[] {childList.getAtom(2),childList.getAtom(3),
                childList.getAtom(4),childList.getAtom(5),childList.getAtom(6),
                childList.getAtom(7),childList.getAtom(8),childList.getAtom(9)};
        dn = childList.getAtom(4);
        dnNon = new IAtom[] {childList.getAtom(3),childList.getAtom(2),childList.getAtom(1),
                childList.getAtom(0)};
        dnLast = childList.getAtom(8);
        dnLastNon = new IAtom[] {childList.getAtom(7),childList.getAtom(6),
                childList.getAtom(5),childList.getAtom(4),childList.getAtom(3),
                childList.getAtom(2),childList.getAtom(1),childList.getAtom(0)};
        iterate = target;
        iterateFirst = targetFirst;
        iterateLast = targetLast;
    }

    /**
     * Test AdjacentPairIterator.
     * @param api the iterator
     */
    private void adjacentPairTests(ApiIntragroup api) {

        api.setBasis(new MoleculeSetSinglet(parent));

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiTwoIterates(api, new AtomPair(iterate, up), new AtomPair(dn, iterate));

        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testApiOneIterate(api, new AtomPair(iterateFirst, upFirst));

        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testApiOneIterate(api, new AtomPair(dnLast, iterateLast));

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testApiOneIterate(api, new AtomPair(dnLast, iterateLast));

        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testApiOneIterate(api, new AtomPair(iterateFirst, upFirst));


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testApiOneIterate(api, new AtomPair(iterate, up));

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testApiOneIterate(api, new AtomPair(dn, iterate));

        //no target, n-1 iterates
        api.setTarget(null);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), parent.getChildList().getAtomCount()-1);

        //if no target, direction doesn't matter
        api.setDirection(null);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        //no basis, no iterates
        api.setBasis(null);
        testNoIterates(api);
        api.setTarget(target);
        testNoIterates(api);
    }

    /**
     * Test NonAdjacentPairIterator.
     * @param api the iterator
     */
    private void nonAdjacentPairTests(ApiIntragroup api) {

        api.setBasis(new MoleculeSetSinglet(parent));

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiIterates(api, iterate, upNon, dnNon);

        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testApiIterates(api, UP, iterateFirst, upFirstNon);

        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testApiIterates(api, DOWN, iterateLast, dnLastNon);

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testApiIterates(api, DOWN, iterateLast, dnLastNon);

        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testApiIterates(api, UP, iterateFirst, upFirstNon);


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testApiIterates(api, UP, iterate, upNon);

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testApiIterates(api, DOWN, iterate, dnNon);

        //no target, (2(n-2) + (n-2)*(n-3))/2 iterates
        //(first and last have n-2, the other n-2 have n-3)
        api.setTarget(null);
        LinkedList list0 = generalIteratorMethodTests(api);
        int n = parent.getChildList().getAtomCount();
        assertEquals(list0.size(), (2*(n-2) + (n-2)*(n-3))/2);

        //if no target, direction doesn't matter
        api.setDirection(null);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        //no basis, no iterates
        api.setBasis(null);
        testNoIterates(api);
        api.setTarget(target);
        testNoIterates(api);
    }
}
