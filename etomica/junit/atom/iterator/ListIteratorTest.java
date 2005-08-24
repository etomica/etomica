package etomica.junit.atom.iterator;
import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomListTabbed;
import etomica.junit.UnitTest;
import etomica.space3d.Space3D;

/**
 * General unit test for a list iterator.  Sets up different types of lists,
 * with different numbers of atoms and tabs, which are then tested for
 * different iterator conditions by subclasses written for particular 
 * list iterators.
 *  
 * @author Ken Benjamin
 * 
 */

/*
 * Created on Oct 1, 2004
 */

public abstract class ListIteratorTest extends IteratorTest {
	
	public ListIteratorTest() {
        tabType1 = AtomLinker.Tab.requestTabType();
        tabType2 = AtomLinker.Tab.requestTabType();
	}

    /**
     * Tests iterator in the various conditions it may have, such as its direction,
     * starting/ending element, and other such features as they are available
     * in the iterator being tested.  This is invoked by testListVariations for
     * different types of lists.
     */
	public abstract void iteratorStateTests(AtomList list);

    /**
     * Sets up different lists, changing the number of atoms and tabs,
     * and for each invokes iteratorStateTests for checking.
     */
	public void testListVariations() {
		
		Space space = Space3D.getInstance();
//		 make empty list to start
		AtomList atomList = new AtomListTabbed();
		
		iteratorStateTests(atomList); 

        AtomLinker.Tab tab1 = AtomLinker.newTab(atomList, tabType1);
        AtomLinker.Tab tab2 = AtomLinker.newTab(atomList, tabType2);
        AtomLinker.Tab tab3 = AtomLinker.newTab(atomList, tabType1);
		for(int i=1; i<11; i++) {
			atomList.add(new Atom(space));
            //test with no tabs
            iteratorStateTests(atomList); 
            for(int j=0; j<i; j++) {
                //test with one tab, testing at every position in the list
                atomList.addBefore(tab1, atomList.entry(j));
                iteratorStateTests(atomList);
                for(int k=0; k<i; k++) {
                    //test with two different tabs, each at every position in list
                    atomList.addBefore(tab2, atomList.entry(k));
                    //test without considering 2nd tab
                    iteratorStateTests(atomList);
                    //test while considering 2nd tab
                    iteratorStateTests(atomList);
                    for(int m=0; m<i; m++) {
                        //test with three different tabs, each at every position in list
                        atomList.addBefore(tab3, atomList.entry(m));
//                        System.out.println(atomList.size()+" "+i+" "+j+" "+k+" "+m);
                        iteratorStateTests(atomList);
                        atomList.remove(tab3);
                    }
                    atomList.remove(tab2);
                }
                atomList.remove(tab1);
            }
        }
        
        //		put a tab at the beginning of the list
		atomList.addFirst(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(atomList);

        //		put a tab at the end of the list		
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(atomList);

//		put a tab in the middle of the list
		AtomLinker.Tab newTab = AtomLinker.newTab(atomList, tabType1);
		atomList.addBefore(newTab, atomList.entry(5));
		iteratorStateTests(atomList);

//		put adjacent tabs in the list; second tab of a different type
		newTab = AtomLinker.newTab(atomList, tabType2);
		atomList.addBefore(newTab, atomList.entry(5));
		iteratorStateTests(atomList); 
		if(UnitTest.VERBOSE) System.out.println("The size of the list: "+ atomList.size());

        //set up lists containing just tabs
		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(atomList);
		
		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(atomList);
		
		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType2));
		iteratorStateTests(atomList);
		
	}
    
    protected int tabType1, tabType2;

}
