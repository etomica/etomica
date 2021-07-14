/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.DipoleSourceMolecular;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Dipole;

/**
 * Computes the squared dipole moment in a system.
 */
public class MeterDipoleSumSquared extends DataSourceScalar {

    private final Box box;
    private final DipoleSourceMolecular dipoleSource;

	public MeterDipoleSumSquared(Box box, DipoleSourceMolecular dipoleSource) {
		super("TIP4P water, dipoleSum^2", new CompoundDimension(new Dimension[]{Dipole.DIMENSION},new double[]{2.0}));
		this.dipoleSource = dipoleSource;
		this.box=box;
	}
	public double getDataAsScalar() {
		Vector dipoleSum = box.getSpace().makeVector();
		IMoleculeList moleculeList = box.getMoleculeList();
		int numMolecule = moleculeList.size();
		for (int i=0;i<numMolecule; i++){
			dipoleSum.PE(dipoleSource.getDipole(moleculeList.get(i)));
		}
        return dipoleSum.squared();
	}
}
