package etomica.junit;

import etomica.*;
import junit.framework.*;

/**
 * @author aawalker
 *
 */
abstract class IteratorTest extends TestCase {
	
	protected Simulation sim;
	protected Lister[] lister;
	protected Species species1,species2,species3;
	protected int nMolecules=15;
	protected int nAtoms=3;
	protected int nLists=8; //minimum 8
	protected Phase phase;
	protected Atom basis1,basis2;
	protected AtomPair pairBasis;
	protected AtomList atomList,atomList1,atomList2;
	protected AtomLinker.Tab header;
	protected AtomLinker firstLinker1,lastLinker1,midLinker1,firstLinker2,lastLinker2,midLinker2;
	protected Atom first1,middle1,last1,first2,middle2,last2,outsideAtom;
	protected int mid1,mid2;
	
	/**
	 * This method is called before any test is run and sets up the basics of the
	 * etomica simulation and testing elements.
	 */
	protected void setUp() {
		sim = new Simulation();
		phase = new Phase();
		species1 = new SpeciesSpheresMono(nMolecules);
		species2 = new SpeciesSpheresMono(nMolecules);
		species3 = new SpeciesSpheresMono(1);
		sim.elementCoordinator.go();
		lister = Lister.listerArray(nLists,sim.space());
		basis1 = species1.getAgent(phase);
		basis2 = species2.getAgent(phase);
		pairBasis = new AtomPair(basis1,basis2);
		atomList = phase.speciesMaster.atomList;
		atomList1 = ((AtomTreeNodeGroup)basis1.node).childList;
		atomList2 = ((AtomTreeNodeGroup)basis2.node).childList;
		first1 = atomList1.getFirst();
		last1 = atomList1.getLast();
		mid1=atomList1.size()/2;
		middle1 = atomList1.get(mid1);
		first2 = atomList2.getFirst();
		last2 = atomList2.getLast();
		mid2 = atomList2.size()/2;
		middle2 = atomList2.get(mid2);
		header = atomList.header;
		firstLinker1 = atomList1.firstEntry();
		lastLinker1 = atomList1.lastEntry();
		midLinker1 = atomList1.entry(mid1);
		firstLinker2 = atomList2.firstEntry();
		lastLinker2 = atomList2.lastEntry();
		midLinker2 = atomList2.entry(mid2);
		outsideAtom = species3.getAgent(phase).firstMolecule();
	}

	/**
	 * Clears all of the lists.
	 */
		protected void clearLists() {
			for (int i=0;i<nLists;i++) {
				lister[i].list.clear();
			}
		}
	
	/**
	 * Prints all of the lists.
	 */
		protected void printLists() {
			for (int i=0;i<nLists;i++) {
				System.out.println(lister[i].list);
			}
			System.out.println();
		}
}
