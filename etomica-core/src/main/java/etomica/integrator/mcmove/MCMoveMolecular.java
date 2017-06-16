/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.box.Box;
import etomica.molecule.iterator.MoleculeIterator;



/**
 * 
 * An interface that contains MoleculeIterator
 * 
 * @author taitan
 *
 */
public interface MCMoveMolecular{

	public MoleculeIterator affectedMolecules(Box box);
	
}
