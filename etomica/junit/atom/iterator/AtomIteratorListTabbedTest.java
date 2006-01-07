package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomListTabbed;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.space3d.Space3D;

/**
 * Unit test for AtomIteratorList. Tests that effect of first and terminator type are
 * correct, and checks for appropriate behavior with or without tabs in list.
 * 
 * @author Ken Benjamin, David Kofke
 *  
 */

/*
 * Created on Oct 1, 2004
 */

public class AtomIteratorListTabbedTest extends IteratorTestAbstract {

    protected final int nLists = 10;
    protected LinkedList[] lists = new LinkedList[nLists];

    public AtomIteratorListTabbedTest() {
        super();
    }
    
    /**
     * Test that documented exceptions are thrown.
     */
    public void testExceptions() {
        boolean exceptionThrown = false;
        Space space = Space3D.getInstance();

        int tabType1 = AtomLinker.Tab.requestTabType();
        int tabType2 = AtomLinker.Tab.requestTabType();
       
        AtomList atomList = new AtomListTabbed();
        
        AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(atomList);
        AtomLinker.Tab tab1 = AtomLinker.newTab(atomList, tabType1);
        AtomLinker.Tab tab2 = AtomLinker.newTab(new AtomList(), tabType1);
        AtomLinker.Tab tab3 = AtomLinker.newTab(atomList, tabType2);
        
        for(int i=0; i<3; i++) {
            atomList.add(new Atom(space));
        }
        atomList.add(tab1);
        
        try {
            iterator.setFirst(tab3);
        } catch(IllegalArgumentException ex) {
            exceptionThrown = true;
        }
        assertFalse(exceptionThrown);

        
        try {
            iterator.setFirst(tab2);
        } catch(IllegalArgumentException ex) {
            exceptionThrown = true;
        }
        
        assertTrue(exceptionThrown);
    }

    /**
     * Sets up different lists, changing the number of atoms and tabs,
     * and for each invokes iteratorStateTests for checking.
     */
    public void testListVariations() {
        
        Space space = Space3D.getInstance();
        int tabType1 = AtomLinker.Tab.requestTabType();
        int tabType2 = AtomLinker.Tab.requestTabType();

        
        // *** ZERO TABS ***
        AtomList atomList = new AtomList();

        //test empty list
        tab0Tests(atomList);
        
        //test list with atoms but no tabs
        for(int i=1; i<11; i++) {
            atomList.add(new Atom(space));
            tab0Tests(atomList);
        }
        
        // *** ONE TAB ***
        atomList = new AtomListTabbed();
        AtomLinker.Tab tab1 = AtomLinker.newTab(atomList, tabType1);
        
        //test empty list with one tab
        atomList.clear();
        atomList.add(tab1);
        tab1Tests(atomList, tab1, 0);
        atomList.remove(tab1);
        
        //test list with atoms and one tab
        for(int i=1; i<5; i++) {
            atomList.add(new Atom(space));
            for(int j=0; j<i; j++) {
                atomList.addBefore(tab1, atomList.entry(j));
                tab1Tests(atomList, tab1, j);
                atomList.remove(tab1);
            }
            atomList.add(tab1);//add last
            tab1Tests(atomList, tab1, atomList.size());
            atomList.remove(tab1);
        }
        
        // *** TWO TABS, SAME TYPE ***
        AtomLinker.Tab tab2 = AtomLinker.newTab(atomList, tabType1);
        tab2Tests(tab1, tab2, tab2SametypeTester);

        // *** TWO TABS, DIFFERENT TYPE ***
        tab2 = AtomLinker.newTab(atomList, tabType2);
        tab2Tests(tab1, tab2, tab2DifferenttypeTester);        
    }
    
    //generates variation of lists with two tabs and different numbers of atoms
    public void tab2Tests(AtomLinker.Tab tab1, AtomLinker.Tab tab2, Tab2Tester tester) {
        // test empty list with two tabs
        Space space = Space1D.getInstance();
        AtomList atomList = new AtomListTabbed();
        atomList.add(tab1);
        atomList.add(tab2);
        tester.tab2Tests(atomList, tab1, 0, tab2, 0);
        atomList.remove(tab1);
        atomList.remove(tab2);

        //test list with atoms and two tabs
        for(int i=1; i<11; i++) {
            atomList.add(new Atom(space));
            for(int j=0; j<i; j++) {
                atomList.addBefore(tab1, atomList.entry(j));
                for(int k=0; k<i; k++) {
                    atomList.addBefore(tab2, atomList.entry(k));
                    tester.tab2Tests(atomList, tab1, j, tab2, k);
                    atomList.remove(tab2);
                }
                atomList.add(tab2);//tab2 after all atoms
                tester.tab2Tests(atomList, tab1, j, tab2, atomList.size());
                atomList.remove(tab2);
                
                atomList.remove(tab1);
            }
            atomList.add(tab1);//tab1 after all atoms
            for(int k=0; k<i; k++) {
                atomList.addBefore(tab2, atomList.entry(k));
                tester.tab2Tests(atomList, tab1, atomList.size(), tab2, k);
                atomList.remove(tab2);
            }
            atomList.add(tab2);//tab2 after all atoms and tab1
            tester.tab2Tests(atomList, tab1, atomList.size(), tab2, atomList.size());
            atomList.remove(tab2);
            
            atomList.remove(tab1);
            
        }
    }

    private final Tab2Tester tab2SametypeTester = new Tab2Tester() {
        public void tab2Tests(AtomList list, AtomLinker.Tab tab1, int index1,
                AtomLinker.Tab tab2, int index2) {
            if(tab1.type != tab2.type) throw new RuntimeException("test set up improperly");
            
            AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(list);
            iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);
    
            //iterate entire list
            lists[0] = generalIteratorMethodTests(iterator);
            assertEquals(lists[0].size(), list.size());
            if (!list.isEmpty()) {
                checkFirstLast(lists[0], list.getFirst(), list.getLast());
            }
    
            //test iteration starting from tab1, ignoring tab2
            iterator.setFirst(tab1);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), list.size()-index1);
            //test iteration starting from tab2, ignoring tab1
            iterator.setFirst(tab2);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), list.size()-index2);
            
            //test iteration starting from tab1 and ending with their type
            iterator.setFirst(tab1);
            iterator.setTerminatorType(tab1.type);
            lists[1] = generalIteratorMethodTests(iterator);
            int nIter = (index2 >=index1) ? (index2-index1) : list.size()-index1;
            assertEquals(lists[1].size(), nIter);
            iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1], lists[2]);

            //test iteration starting from beginning and ending with their type
            iterator.setFirst(null);
            iterator.setTerminatorType(tab1.type);
            lists[1] = generalIteratorMethodTests(iterator);
            nIter = Math.min(index2,index1);
            assertEquals(lists[1].size(), nIter);
            iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1], lists[2]);
        }
    };//end of tab2SametypeTester

    private final Tab2Tester tab2DifferenttypeTester = new Tab2Tester() {
        public void tab2Tests(AtomList list, AtomLinker.Tab tab1, int index1,
                AtomLinker.Tab tab2, int index2) {
            if(tab1.type == tab2.type) throw new RuntimeException("test set up improperly");
            
            AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(list);
            iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);
    
            //iterate entire list
            lists[0] = generalIteratorMethodTests(iterator);
            assertEquals(lists[0].size(), list.size());
            if (!list.isEmpty()) {
                checkFirstLast(lists[0], list.getFirst(), list.getLast());
            }
    
            //test iteration starting from tab1, ignoring tab2
            iterator.setFirst(tab1);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), list.size()-index1);
            //test iteration starting from tab2, ignoring tab1
            iterator.setFirst(tab2);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), list.size()-index2);
            
            //test iteration starting from tab1 and ending with tab1 type
            iterator.setFirst(tab1);
            iterator.setTerminatorType(tab1.type);
            lists[1] = generalIteratorMethodTests(iterator);
            int nIter = list.size()-index1;
            assertEquals(lists[1].size(), nIter);
            //test iteration starting from tab2 and ending with tab2 type
            iterator.setFirst(tab2);
            iterator.setTerminatorType(tab2.type);
            lists[1] = generalIteratorMethodTests(iterator);
            nIter = list.size()-index2;
            assertEquals(lists[1].size(), nIter);
            //test iteration starting from tab1 and ending with tab2 type
            iterator.setFirst(tab1);
            iterator.setTerminatorType(tab2.type);
            lists[1] = generalIteratorMethodTests(iterator);
            nIter = (index2 >=index1) ? (index2-index1) : list.size()-index1;
            assertEquals(lists[1].size(), nIter);
            iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1], lists[2]);
            //test iteration starting from tab2 and ending with tab1 type
            iterator.setFirst(tab2);
            iterator.setTerminatorType(tab1.type);
            lists[1] = generalIteratorMethodTests(iterator);
            //note difference from reverse -- tab2 is after tab1 in list when adjacent
            nIter = (index1 >index2) ? (index1-index2) : list.size()-index2;
            assertEquals(lists[1].size(), nIter);
            iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1], lists[2]);
 
            //test iteration starting from beginning and ending with tab1 type
            iterator.setFirst(null);
            iterator.setTerminatorType(tab1.type);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), index1);

            //test iteration starting from beginning and ending with tab2 type
            iterator.setTerminatorType(tab2.type);
            lists[1] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1].size(), index2);

            iterator.setFirst(null);
            iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);
            lists[1] = generalIteratorMethodTests(iterator);
            nIter = Math.min(index2,index1);
            assertEquals(lists[1].size(), nIter);
            
            //test that zero gives header as terminator
            iterator.setTerminatorType(0);
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[0], lists[2]);

            //test that negative value gives any tab as terminator
            iterator.setTerminatorType(-100);
            lists[2] = generalIteratorMethodTests(iterator);
            assertEquals(lists[1], lists[2]);

        }
    };//end of tab2DifferenttypeTester

    /**
     * Performs tests on a list having one tab
     * @param list the list
     * @param tab the tab
     * @param index the index of the atom following the tab
     */
    public void tab1Tests(AtomList list, AtomLinker.Tab tab, int index) {
        AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(list);
        iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);

        //iterate entire list
        lists[0] = generalIteratorMethodTests(iterator);
        assertEquals(lists[0].size(), list.size());
        if (!list.isEmpty()) {
            checkFirstLast(lists[0], list.getFirst(), list.getLast());
        }

        //test iteration starting from tab
        iterator.setFirst(tab);
        lists[1] = generalIteratorMethodTests(iterator);
        assertEquals(lists[1].size(), list.size()-index);
        iterator.setTerminatorType(tab.type);//shouldn't change anything
        lists[2] = generalIteratorMethodTests(iterator);
        assertEquals(lists[1], lists[2]);
        iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
        lists[2] = generalIteratorMethodTests(iterator);
        assertEquals(lists[1], lists[2]);

        //test whether can set iteration over whole list again
        iterator.setFirst(null);
        iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);
        lists[1] = generalIteratorMethodTests(iterator);
        assertEquals(lists[0], lists[1]);
        
        //test iteration ending at tab
        iterator.setTerminatorType(tab.type);
        lists[1] = generalIteratorMethodTests(iterator);
        assertEquals(lists[1].size(), index);
        iterator.setTerminatorType(AtomLinker.Tab.ANY_TAB);//shouldn't change anything
        lists[2] = generalIteratorMethodTests(iterator);
        assertEquals(lists[1], lists[2]);
        
        //test that zero terminatorType corresponds to header
        iterator.setTerminatorType(0);
        lists[1] = generalIteratorMethodTests(iterator);
        assertEquals(lists[0], lists[1]);
    }

    /**
     * Performs test on a list having no tabs
     * @param list the list
     */
    public void tab0Tests(AtomList list) {
        AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(list);
        iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);

        //iterate entire list
        lists[0] = generalIteratorMethodTests(iterator);
        assertEquals(lists[0].size(), list.size());
        if (!list.isEmpty()) {
            checkFirstLast(lists[0], list.getFirst(), list.getLast());
        }
    }

    private void checkFirstLast(LinkedList list, Atom first, Atom last) {
        if (list.size() > 0) {
            assertEquals(list.getFirst(), first.toString());
            assertEquals(list.getLast(), last.toString());
        }
    }

    //performs tests on iterator with two tabs.  Instances are defined
    //here to handle cases in which tabs are the same or different type
    private interface Tab2Tester {
        public void tab2Tests(AtomList list, AtomLinker.Tab tab1, int index1,
                AtomLinker.Tab tab2, int index2);
    }
}

