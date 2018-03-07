/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

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
}
