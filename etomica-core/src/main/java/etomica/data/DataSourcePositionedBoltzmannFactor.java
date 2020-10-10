/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Null;

/**
 * Calculates the Boltzmann factor at a position within a box a molecule of a
 * particular species would have if it existed at that point.  Only works for
 * monatomic molecules (no rotation is attempted).
 * 
 * @author Andrew Schultz
 */
public class DataSourcePositionedBoltzmannFactor implements DataSourcePositioned {

    public DataSourcePositionedBoltzmannFactor(Space space, SpeciesGeneral species) {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("chemical potential", Null.DIMENSION);
        tag = new DataTag();
        atomTranslator = new MoleculeActionTranslateTo(space);
        this.species = species;
    }

    /**
     * Sets the integrator.  The integrator is used to obtain the
     * PotentialMaster, Box and temperature.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
        energyMeter = new MeterPotentialEnergy(integrator.getPotentialMaster(), integrator.getBox());
    }

    /**
     * Returns the ISpecies for which the chemical potential is to be measured.
     */
    public ISpecies getSpecies() {
        return species;
    }

    public IData getData(Vector a) {
        Box box = integrator.getBox();
        atomTranslator.setDestination(a);
        IMolecule testMolecule = box.addNewMolecule(species, mol -> {
            atomTranslator.actionPerformed(mol);
        });
        energyMeter.setTarget(testMolecule);
        double temp = integrator.getTemperature();
        data.x = Math.exp(-energyMeter.getDataAsScalar()/temp);

        box.removeMolecule(testMolecule);

        return data;
    }

    public void setBox(Box newBox) {
        if (integrator != null && newBox != integrator.getBox()) {
            throw new RuntimeException("You should really figure out which Box is for me!");
        }
    }

    public IDataInfo getPositionDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected MoleculeActionTranslateTo atomTranslator;
    protected MeterPotentialEnergy energyMeter;
    protected final DataInfoDouble dataInfo;
    protected final DataDouble data;
    protected IntegratorBox integrator;
    protected final DataTag tag;
    private final ISpecies species;
}
