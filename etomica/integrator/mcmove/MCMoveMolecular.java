package etomica.integrator.mcmove;

import etomica.api.IBox;
import etomica.atom.iterator.MoleculeIterator;



/**
 * 
 * An interface that contains MoleculeIterator
 * 
 * @author taitan
 *
 */
public interface MCMoveMolecular{

	public MoleculeIterator affectedMolecules(IBox box);
	
}
