package etomica.junit;

import java.util.LinkedList;

import etomica.Atom;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorList;

/**
 * @author aawalker
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class ApiIntergroup1ATest extends IteratorTest {
	
	private ApiIntergroup1A iterator;
	
	public static void main(String[] args) {
		junit.textui.TestRunner.run(ApiIntergroup1A.class);
	}
	
/**
 * This method is called before any test is run and sets up the basics of the
 * etomica simulation and testing elements.
 */
	protected void setUp() {
		super.setUp();
		iterator = new ApiIntergroup1A(sim);
	}
		

/**
 * returns the list expected from an ApiIntergroup1A iteration.
 * 
 * @param singleAtom the single atom which is a member of all pairs
 * @param group the group of atoms which is being paired with singleAtom
 * @return LinkedList a list of the AtomPairs toString method
 */	
	private LinkedList lister1A(Atom singleAtom, Atom group) {
		LinkedList list = new LinkedList();
		AtomPair ap;
		AtomIteratorList listIterator = new AtomIteratorList(((AtomTreeNodeGroup)group.node).childList);
		while (listIterator.hasNext()) {
			list.add((new AtomPair(singleAtom,listIterator.next())).toString());
		}
		return list;
	}

/**
 * Tests the basic has/next iteration
 */	
	public void testHasNextIteration() {

		iterator.setBasis(basis1,basis2);

		iterator.reset(first1);		
		while (iterator.hasNext()) {
			lister[0].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list,lister1A(first1,basis2));
		
		iterator.reset(first2);
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(lister[1].list,lister1A(first2,basis1));
//		System.out.println(lister[0].list);
//		System.out.println(lister[1].list);

		iterator.reset(outsideAtom);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(lister[2].list,new LinkedList());
		
		clearLists();
	}
	
}
