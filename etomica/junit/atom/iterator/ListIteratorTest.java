package etomica.junit.atom.iterator;
import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
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
		
		Space space = new Space3D();
//		 make empty list to start
		AtomList atomList = new AtomList();
		
//		listElements();
		int nMolecules=0;
		if(UnitTest.VERBOSE) System.out.println("The size of the empty list is "+atomList.size());
		if(UnitTest.VERBOSE) System.out.println("The element in the empty list is "+atomList.getFirst());
		iteratorStateTests(atomList); 

//		set zeroth element of array to null
/*		atomList.add(0, null);
		listElements();
		iterator.setList(atomList);
		iteratorStateTests(iterator);
		System.out.println("Just created list with null element in middle");
*/// 		added first atom; 1 element in list
		atomList.add(new Atom(space));
		nMolecules=++nMolecules;
		if(UnitTest.VERBOSE) System.out.println("The value of nMolecules is: "+nMolecules);
		iteratorStateTests(atomList);  
//		 added second atom; 2 elements in list
		atomList.add(new Atom(space));
		nMolecules=++nMolecules;
		if(UnitTest.VERBOSE) System.out.println("The value of nMolecules is: "+nMolecules);
		iteratorStateTests(atomList);
//		 adding eight atoms; 10 elements in list
		for(int i=0; i<8; i++) {
			atomList.add(new Atom(space));
			nMolecules=++nMolecules;
			if(UnitTest.VERBOSE) System.out.println("The value of nMolecules is: "+nMolecules);
		}
		iteratorStateTests(atomList);
		if(UnitTest.VERBOSE) System.out.println("The size of the list: "+ atomList.size());

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
