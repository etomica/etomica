/*
 * Created on Dec 3, 2004
 */
package etomica.models.hexane;

import etomica.config.ConfigurationLattice;
import etomica.space.Space;


/**
 * @author nancycribbin
 *
 * Creates the lattice that the hexane molecules of Dr. Monson's data are placed
 * into.
 */
public class ConfigurationHexane extends ConfigurationLattice {

	ConfigurationHexane(Space space) {
	    super(new PrimitiveHexane(space));
	}
}