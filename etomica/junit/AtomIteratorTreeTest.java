package etomica.junit;

import etomica.*;
import etomica.action.AtomAction;
import junit.framework.*;
import java.util.*;

/**
 * @author aawalker
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AtomIteratorTreeTest extends TestCase {

	private Species species0,species1;
	private AtomIteratorTree iterator;
	private AtomIteratorList listIterator;
	private IteratorDirective id;
	private final int nMolecules=10;
	private final int nAtoms0=2;
	private final int nAtoms1=3;
	private final int nLists=8;
	private Phase phase;
	private LinkedList[] list;
	private AtomAction[] lister;
	private Atom basis0,basis1,basis2;
	private AtomList atomList0,atomList1,atomList2;

	public static void main(String[] args) {
		junit.textui.TestRunner.run(AtomIteratorTreeTest.class);
	}

	protected void setUp() {
		Simulation sim = new Simulation();
		phase = new Phase();
		species0 = new SpeciesSpheres(nMolecules/2,2);
		species1 = new SpeciesSpheres(nMolecules/2,3);
		sim.elementCoordinator.go();
		basis0 = species0.getAgent(phase);
		basis1 = species1.getAgent(phase);
		basis2 = phase.speciesMaster;
		atomList0 = ((AtomTreeNodeGroup)species0.getAgent(phase).node).childList;
		atomList1 = ((AtomTreeNodeGroup)species1.getAgent(phase).node).childList;
		atomList2 = phase.speciesMaster.atomList;
		iterator = new AtomIteratorTree();
		listIterator = new AtomIteratorList();
		id = new IteratorDirective();
		list=new LinkedList[nLists];
		for (int i=0;i<nLists;i++) {
			list[i]=new LinkedList();
		}
		lister=new AtomAction[nLists];
		for (int i=0;i<list.length;i++) {
			lister[i]=new Lister(i);
			}
		}
	
	private void clearLists() {
		for (int i=0;i<nLists;i++) {
			list[i].clear();
		}
	}
	
	private void printLists() {
		for (int i=0;i<nLists;i++) {
			System.out.println(list[i]);
		}
		System.out.println();
	}

	public void testNMomlecules() {

		iterator.setRoot(basis2);

		iterator.setIterationDepth(0);

		iterator.setRoot(basis0);

		assertEquals(1,iterator.size());
		
		iterator.setRoot(basis1);
		assertEquals(1,iterator.size());
		
		iterator.setRoot(basis2);
		assertEquals(1,iterator.size());
		
		iterator.setIterationDepth(1);
		
		iterator.setRoot(basis0);
		assertEquals(nMolecules/2,iterator.size());
//		System.out.println(iterator.size());

		iterator.setRoot(basis1);
		assertEquals(nMolecules/2,iterator.size());
//		System.out.println(iterator.size());

		iterator.setRoot(basis2);
		assertEquals(2,iterator.size());
//		System.out.println(iterator.size());
		
		iterator.setIterationDepth(2);

		iterator.setRoot(basis0);
		assertEquals((nMolecules/2)*nAtoms0,iterator.size());
		
		iterator.setRoot(basis1);
		assertEquals((nMolecules/2)*nAtoms1,iterator.size());
		
		iterator.setRoot(basis2);
		assertEquals(nMolecules,iterator.size());
		
		iterator.setIterationDepth(3);

		iterator.setRoot(basis0);
		assertEquals((nMolecules/2)*nAtoms0,iterator.size());
		
		iterator.setRoot(basis1);
		assertEquals((nMolecules/2)*nAtoms1,iterator.size());
		
		iterator.setRoot(basis2);
		assertEquals((nMolecules/2)*(nAtoms0+nAtoms1),iterator.size());

		iterator.setAsLeafIterator();
		
		iterator.setRoot(basis0);
		assertEquals((nMolecules/2)*nAtoms0,iterator.size());
		
		iterator.setRoot(basis1);
		assertEquals((nMolecules/2)*nAtoms1,iterator.size());
		
		iterator.setRoot(basis2);
		assertEquals((nMolecules/2)*(nAtoms0+nAtoms1),iterator.size());		
	}
	
	public void testAllDepth0() {

		iterator.setIterationDepth(0);
		
		iterator.all(basis0,id,lister[0]);
		LinkedList list0 = new LinkedList();
		list0.add(basis0.toString());
		assertEquals(list0,list[0]);
		
		clearLists();
	}
	
	public void testAllDepth1() {
		
		iterator.setIterationDepth(1);

		iterator.all(basis0,id,lister[0]);
		listIterator.all(basis0,id,lister[1]);
		assertEquals(list[0],list[1]);
		
		iterator.all(basis1,id,lister[2]);
		listIterator.all(basis1,id,lister[3]);
		assertEquals(list[2],list[3]);
		
		iterator.all(basis2,id,lister[4]);
		listIterator.all(basis2,id,lister[5]);
		assertEquals(list[4],list[5]);
		
		clearLists();
	}
	
	public void testAllLeaf() {
		
		AtomIteratorList listIterator0 = new AtomIteratorList(atomList0);
		AtomIteratorList listIterator1 = new AtomIteratorList(atomList1);

		iterator.setAsLeafIterator();

		iterator.all(basis0,id,lister[0]);
		while (listIterator0.hasNext()) {
			listIterator.all(listIterator0.next(),id,lister[1]);
		}
		assertEquals(list[0],list[1]);
		
		iterator.all(basis2,id,lister[2]);
		while (listIterator1.hasNext()) {
			listIterator.all(listIterator1.next(),id,lister[1]);
		}
//		System.out.println(list[2]);
//		System.out.println(list[1]);
		assertEquals(list[2],list[1]);
		
		clearLists();
	}
	
	


	class Lister implements AtomAction {

		private final int i;

		public Lister(int i) {
			this.i=i;
		}
		public void actionPerformed(AtomSet as) {
			list[i].add(as.toString());	
		}
		public void actionPerformed(Atom a) {
			list[i].add(a.toString());
		}
	}

}
