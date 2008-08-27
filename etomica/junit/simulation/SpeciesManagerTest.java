package etomica.junit.simulation;

import junit.framework.TestCase;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.ISpeciesManager;
import etomica.atom.AtomTypeLeaf;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresHetero;
import etomica.species.SpeciesSpheresMono;

public class SpeciesManagerTest extends TestCase {

	private ISpeciesManager sm;
	ISimulation simulation;
	ISpace space;
	Element element;
	
	public void setUp() {
		space = Space.getInstance(3);
		simulation = new Simulation(space);
		sm = simulation.getSpeciesManager();
		element = new ElementSimple(simulation);
		if(sm == null) {
			System.out.println("species manager is null");
		}
	}

	/*
	 * testAddSpecies
	 */
	public void testAddSpecies() {
		final int numSpecies = 5;
		ISpecies species[] = new ISpecies[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresMono(simulation, space, element);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());
		
		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertSame(species[i], sm.getSpecies(i));
		}
	}

	/*
	 * testRemoveSpecies
	 */
	public void testRemoveSpecies() {
		int numSpecies = 5;
		ISpecies species[] = new ISpecies[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresMono(simulation, space, element);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertSame(species[i], sm.getSpecies(i));
		}
		
		sm.removeSpecies(species[numSpecies-1]);
		numSpecies--;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertSame(species[i], sm.getSpecies(i));
		}

		sm.removeSpecies(species[0]);
		numSpecies--;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 1; i < sm.getSpeciesCount(); i++) {
			assertSame(species[i], sm.getSpecies(i-1));
		}
	}

	/*
	 * testAddSpeciesChildIndex
	 */
	public void testAddSpeciesChildIndex() {
		int numSpecies = 5;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < numSpecies; i++) {
			assertEquals(i, sm.getSpecies(i).getIndex());
		}

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}

		SpeciesSpheresMono newSpecies = new SpeciesSpheresMono(simulation, space, element);
		for(int j = 0; j < numAtomsPerSpecies-1; j++) {
		    newSpecies.addChildType(new AtomTypeLeaf(element));
		}
		sm.addSpecies(newSpecies);
		numSpecies++;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}

	}

	/*
	 * testRemoveFirstSpeciesChildIndex
	 */
	public void testRemoveFirstSpeciesChildIndex() {
		int numSpecies = 6;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < numSpecies; i++) {
			assertEquals(i, sm.getSpecies(i).getIndex());
		}

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}

		sm.removeSpecies(species[0]);
		numSpecies--;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}
	}

	/*
	 * testRemoveSpeciesFromMiddleChildIndex
	 */
	public void testRemoveSpeciesFromMiddleChildIndex() {
		int numSpecies = 6;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < numSpecies; i++) {
			assertEquals(i, sm.getSpecies(i).getIndex());
		}

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}

		sm.removeSpecies(species[2]);
		numSpecies--;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}
	}

	/*
	 * testRemoveLastSpeciesChildIndex
	 */
	public void testRemoveLastSpeciesChildIndex() {
		int numSpecies = 6;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    sm.addSpecies(species[i]);
		}

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < numSpecies; i++) {
			assertEquals(i, sm.getSpecies(i).getIndex());
		}

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}

		sm.removeSpecies(species[5]);
		numSpecies--;

		assertEquals(numSpecies, sm.getSpeciesCount());

		for(int i = 0; i < sm.getSpeciesCount(); i++) {
			assertEquals(numAtomsPerSpecies, sm.getSpecies(i).getChildTypeCount());
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(((numSpecies)+(i*numAtomsPerSpecies) + j),
			    		   sm.getSpecies(i).getChildType(j).getIndex());
			}
		}
	}

}
