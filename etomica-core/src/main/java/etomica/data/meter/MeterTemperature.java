/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeOrientedDynamic;
import etomica.species.SpeciesManager;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Temperature;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition.  The class uses a MeterKineticEnergy by default to calculate
 * the kinetic energy, but any DataSourceScalar can be used for this purpose by
 * calling setKineticEnergyMeter.
 * 
 * If the Simulation is not given, the class will assume that all atoms have
 * only translational degrees of freedom.  If the Simulation is given, this
 * class will examine the ISpecies and calculate the actual number of degrees
 * of freedom (more for oriented atoms or molecules).
 * 
 * @author Andrew Schultz
 */
public class MeterTemperature extends DataSourceScalar {

	protected final SpeciesManager sm;
	private final int dim;
    protected Box box;
    protected DataSourceScalar meterKE;

	public MeterTemperature(Box box, int D) {
		this(null, box, D);
	}

	public MeterTemperature(SpeciesManager sm, Box box, int D) {
		super("Temperature", Temperature.DIMENSION);
		dim = D;
		meterKE = new MeterKineticEnergy(box);
		this.sm = sm;
		this.box = box;
	}

    public void setKineticEnergyMeter(DataSourceScalar meterKineticEnergy) {
        meterKE = meterKineticEnergy;
    }

	public double getDataAsScalar() {
		int totalD = box.getLeafList().size() * dim;
		if (sm != null) {
			totalD = 0;
//	        ISpecies[] species = sim.getSpeciesManager().getSpecies();
			for (int i = 0; i < sm.getSpeciesCount(); i++) {
				int nMolecules = box.getNMolecules(sm.getSpecies(i));
				if (nMolecules > 0) {
					IMolecule molecule = box.getMoleculeList(sm.getSpecies(i)).get(0);
					if (molecule instanceof MoleculeOrientedDynamic) {
						if (Double.isInfinite(sm.getSpecies(i).getMass())) {
							continue;
						}
						totalD += 6 * nMolecules;
					} else {
						IAtomList children = molecule.getChildList();
						if (children.size() == 0 ||
								Double.isInfinite(children.get(0).getType().getMass())) {
	                        continue;
	                    }
                        if (children.get(0).getType() instanceof AtomTypeOriented) {
                            // oriented sphere at this point corresponds to cylindrical symmetry
	                        if (dim == 3) {
	                            totalD += 5*nMolecules*children.size();
	                        }
	                        else { // dim = 2
	                            totalD += 3*nMolecules*children.size();
	                        }
	                    }
	                    else {
	                        totalD += dim*nMolecules*children.size();
	                    }
	                }
	            }
	        }
	    }
		return (2. / totalD) * meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Temperature.DIMENSION;
	}
}
