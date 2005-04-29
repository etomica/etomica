/*
 * Created on Oct 1, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.junit;
import junit.framework.TestCase;
import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListDependent;
import etomica.space3d.Space3D;

/**
 * @author Ken Benjamin
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

public abstract class ListIteratorTest extends TestCase {
	
	protected int tabType1, tabType2;
	private final AtomIteratorListDependent iterator;


	public ListIteratorTest(AtomIteratorListDependent iterator) {
		this.iterator = iterator;
		tabType1 = AtomLinker.Tab.requestTabType();
		tabType2 = AtomLinker.Tab.requestTabType();
	}
/**
 * This method is called before any test is run.  It creates the various lists
 * to be tested and calls the iteratorStateTests method for each case.
 */

	public abstract void iteratorStateTests(AtomIteratorListDependent iterator);

	public void testListVariations() {
		
		Space space = new Space3D();
//		 make empty list to start
		AtomList atomList = new AtomList();
		
//		listElements();
		int nMolecules=0;
		iterator.setList(atomList);
		if(UnitTest.VERBOSE) System.out.println("The size of the empty list is "+iterator.size());
		if(UnitTest.VERBOSE) System.out.println("The element in the empty list is "+atomList.getFirst());
		iteratorStateTests(iterator); 

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
		iterator.setList(atomList); // added during kmb3 list iterator testing 1/7/05
		iteratorStateTests(iterator);  
//		 added second atom; 2 elements in list
		atomList.add(new Atom(space));
		nMolecules=++nMolecules;
		if(UnitTest.VERBOSE) System.out.println("The value of nMolecules is: "+nMolecules);
		iteratorStateTests(iterator);
//		 adding eight atoms; 10 elements in list
		for(int i=0; i<8; i++) {
			atomList.add(new Atom(space));
			nMolecules=++nMolecules;
			if(UnitTest.VERBOSE) System.out.println("The value of nMolecules is: "+nMolecules);
		}
		iteratorStateTests(iterator);
		if(UnitTest.VERBOSE) System.out.println("The size of the list: "+ iterator.size());
//		put a tab at the beginning of the list
//		newTab = AtomLinker.newTab(atomList, tabType);
		atomList.addFirst(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(iterator);
//		put a tab at the end of the list		
//		newTab = AtomLinker.newTab(atomList, tabType);
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(iterator);

/*	Working on this block of code below, as of 12/03/04.  Trying to make sure tabs were actually
 * 	placed in the middle of the atom list.  Eventually trying to set two different types of tabs,
 * 	one of which will signal the iteration to stop, the other will permit it to continue.
 * 
 */

//		put a tab in the middle of the list
		AtomLinker.Tab newTab = AtomLinker.newTab(atomList, tabType1);
		atomList.addBefore(newTab, atomList.entry(5));
		iteratorStateTests(iterator);

//		put adjacent tabs in the list; second tab of a different type
		newTab = AtomLinker.newTab(atomList, tabType2);
		atomList.addBefore(newTab, atomList.entry(5));
		iteratorStateTests(iterator);
		if(UnitTest.VERBOSE) System.out.println("The size of the list: "+ iterator.size());

		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(iterator);
		
		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType1));
		iteratorStateTests(iterator);
		
		atomList.clear();
		atomList.add(AtomLinker.newTab(atomList, tabType2));
		iteratorStateTests(iterator);
		
}

}
