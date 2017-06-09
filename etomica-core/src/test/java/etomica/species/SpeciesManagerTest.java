/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.simulation.Simulation;
import etomica.space.Space;
import junit.framework.TestCase;

public class SpeciesManagerTest extends TestCase {

	Simulation simulation;
	Space space;
	IElement element;
	
	public void setUp() {
		space = Space.getInstance(3);
		simulation = new Simulation(space);
		element = new ElementSimple(simulation);
	}

	/*
	 * testAddSpecies
	 */
	public void testAddSpecies() {
		final int numSpecies = 5;
		ISpecies species[] = new ISpecies[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresMono(space, element);
		    simulation.addSpecies(species[i]);
		}

		assertEquals(numSpecies, simulation.getSpeciesCount());
		
		int expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
			assertSame(species[i], simulation.getSpecies(i));
			assertSame(i, species[i].getIndex());
			for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
			    assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
			    expectedChildIndex++;
			}
		}
	}

	/*
	 * testRemoveSpecies
	 */
	public void testRemoveSpecies() {
		int numSpecies = 5;
		ISpecies species[] = new ISpecies[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresMono(space, element);
		    simulation.addSpecies(species[i]);
		}

		assertEquals(numSpecies, simulation.getSpeciesCount());

		int expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
			assertSame(species[i], simulation.getSpecies(i));
			assertSame(i, species[i].getIndex());
	        for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
	            assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
	            expectedChildIndex++;
	        }
		}
		
		simulation.removeSpecies(species[numSpecies-1]);
		numSpecies--;

		assertEquals(numSpecies, simulation.getSpeciesCount());

		expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
			assertSame(species[i], simulation.getSpecies(i));
			assertSame(i, species[i].getIndex());
            for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
                expectedChildIndex++;
            }
		}

		simulation.removeSpecies(species[0]);
		numSpecies--;

		assertEquals(numSpecies, simulation.getSpeciesCount());

		expectedChildIndex = 0;
		for(int i = 1; i < simulation.getSpeciesCount(); i++) {
			assertSame(species[i], simulation.getSpecies(i-1));
			assertSame(i-1, species[i].getIndex());
            for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
                expectedChildIndex++;
            }
		}
	}

	/*
	 * testAddSpeciesChildIndex
	 */
	public void testAddSpeciesChildIndex() {
	    int INIT_NUM_SPECIES = 5;
		int numSpecies = INIT_NUM_SPECIES;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    simulation.addSpecies(species[i]);
		}

		assertEquals(numSpecies, simulation.getSpeciesCount());

		int expectedChildIndex = 0;
		for(int i = 0; i < numSpecies; i++) {
			assertEquals(i, simulation.getSpecies(i).getIndex());
            for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
                expectedChildIndex++;
            }
		}

		SpeciesSpheresMono newSpecies = new SpeciesSpheresMono(space, element);
		for(int j = 0; j < numAtomsPerSpecies-1; j++) {
            newSpecies.addChildType(new AtomType(element));
        }
		simulation.addSpecies(newSpecies);
		numSpecies++;

		assertEquals(numSpecies, simulation.getSpeciesCount());

	    expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
		    if(i < INIT_NUM_SPECIES) {
		        assertEquals(i, species[i].getIndex());
		    }
		    else {
		        assertEquals(i, newSpecies.getIndex());
		    }
			for(int j = 0; j < numAtomsPerSpecies; j++) {
			    assertEquals(expectedChildIndex, simulation.getSpecies(i).getAtomType(j).getIndex());
			    expectedChildIndex++;
			}
		}

	}

	/*
	 * testRemoveFirstSpeciesChildIndex
	 */
	public void testRemoveFirstSpeciesChildIndex() {
        int numSpecies = 6;
        int REMOVE_INDEX = 0;
        final int numAtomsPerSpecies = 3;
        SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

        for(int i = 0; i < numSpecies; i++) {
            species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
            simulation.addSpecies(species[i]);
        }

        assertEquals(numSpecies, simulation.getSpeciesCount());

        int expectedChildIndex = 0;
        for(int i = 0; i < numSpecies; i++) {
            assertEquals(i, species[i].getIndex());
            for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
                expectedChildIndex++;
            }
        }

        simulation.removeSpecies(species[REMOVE_INDEX]);
        numSpecies--;

        assertEquals(numSpecies, simulation.getSpeciesCount());

        expectedChildIndex = 0;
        for(int i = 0; i < simulation.getSpeciesCount(); i++) {
            if(i < REMOVE_INDEX) {
                assertEquals(i, species[i].getIndex());
                for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                    assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
                    expectedChildIndex++;
                }
            }
            else {
                assertEquals(i, species[i+1].getIndex());
                for(int j = 0; j < species[i+1].getAtomTypeCount(); j++) {
                    assertEquals(expectedChildIndex, species[i+1].getAtomType(j).getIndex());
                    expectedChildIndex++;
                }
            }
        }
    }

	/*
	 * testRemoveSpeciesFromMiddleChildIndex
	 */
	public void testRemoveSpeciesFromMiddleChildIndex() {
		int numSpecies = 6;
		int REMOVE_INDEX = 2;
		final int numAtomsPerSpecies = 3;
		SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
		    simulation.addSpecies(species[i]);
		}

		assertEquals(numSpecies, simulation.getSpeciesCount());

		int expectedChildIndex = 0;
		for(int i = 0; i < numSpecies; i++) {
		    assertEquals(i, species[i].getIndex());
			for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
			    assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
			    expectedChildIndex++;
			}
		}

		simulation.removeSpecies(species[REMOVE_INDEX]);
		numSpecies--;

		assertEquals(numSpecies, simulation.getSpeciesCount());

		expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
		    if(i < REMOVE_INDEX) {
		        assertEquals(i, species[i].getIndex());
		        for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
		            assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
		            expectedChildIndex++;
		        }
		    }
		    else {
		        assertEquals(i, species[i+1].getIndex());
	            for(int j = 0; j < species[i+1].getAtomTypeCount(); j++) {
	                assertEquals(expectedChildIndex, species[i+1].getAtomType(j).getIndex());
	                expectedChildIndex++;
	            }
		    }

		}
	}

	/*
	 * testRemoveLastSpeciesChildIndex
	 */
    public void testRemoveLastSpeciesChildIndex() {
        int numSpecies = 6;
        int REMOVE_INDEX = 5;
        final int numAtomsPerSpecies = 3;
        SpeciesSpheresHetero species[] = new SpeciesSpheresHetero[numSpecies];

        for(int i = 0; i < numSpecies; i++) {
            species[i] = new SpeciesSpheresHetero(simulation, space, numAtomsPerSpecies);
            simulation.addSpecies(species[i]);
        }

        assertEquals(numSpecies, simulation.getSpeciesCount());

        int expectedChildIndex = 0;
        for(int i = 0; i < numSpecies; i++) {
            assertEquals(i, species[i].getIndex());
            for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
                expectedChildIndex++;
            }
        }

        simulation.removeSpecies(species[REMOVE_INDEX]);
        numSpecies--;

        assertEquals(numSpecies, simulation.getSpeciesCount());

        expectedChildIndex = 0;
        for(int i = 0; i < simulation.getSpeciesCount(); i++) {
            if(i < REMOVE_INDEX) {
                assertEquals(i, species[i].getIndex());
                for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
                    assertEquals(expectedChildIndex, species[i].getAtomType(j).getIndex());
                    expectedChildIndex++;
                }
            }
            else {
                assertEquals(i, species[i+1].getIndex());
                for(int j = 0; j < species[i+1].getAtomTypeCount(); j++) {
                    assertEquals(expectedChildIndex, species[i+1].getAtomType(j).getIndex());
                    expectedChildIndex++;
                }
            }
        }
    }

}
