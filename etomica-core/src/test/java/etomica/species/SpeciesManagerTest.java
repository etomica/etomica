/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.simulation.Simulation;
import etomica.space.Space;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class SpeciesManagerTest {

	Simulation simulation;
	Space space;
	IElement element;
	
	@BeforeEach
	public void setUp() {
		space = Space.getInstance(3);
		simulation = new Simulation(space);
		element = new ElementSimple(simulation);
	}

	/*
	 * testAddSpecies
	 */
	@Test
	public void testAddSpecies() {
		final int numSpecies = 5;
		ISpecies species[] = new ISpecies[numSpecies];

		for(int i = 0; i < numSpecies; i++) {
		    species[i] = new SpeciesSpheresMono(space, element);
		    simulation.addSpecies(species[i]);
		}

		Assertions.assertEquals(numSpecies, simulation.getSpeciesCount());
		
		int expectedChildIndex = 0;
		for(int i = 0; i < simulation.getSpeciesCount(); i++) {
			Assertions.assertSame(species[i], simulation.getSpecies(i));
			Assertions.assertSame(i, species[i].getIndex());
			for(int j = 0; j < species[i].getAtomTypeCount(); j++) {
			    Assertions.assertSame(expectedChildIndex, species[i].getAtomType(j).getIndex());
			    expectedChildIndex++;
			}
		}
	}
}
