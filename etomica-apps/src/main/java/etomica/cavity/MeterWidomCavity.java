/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

/**
 *
 */
public class MeterWidomCavity implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tag, rTag;
    protected DataDoubleArray rData;
    protected DataDoubleArray.DataInfoDoubleArray rDataInfo;
    protected final IRandom random;
    private int nInsert;
    private ISpecies species;
    private final MeterPotentialEnergy energyMeter;
    protected Box box;
    protected final Vector dr, r0;

    public MeterWidomCavity(Box box, IRandom random, PotentialMaster potentialMaster) {
        setNInsert(100);
        this.box = box;
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster, box);
        tag = new DataTag();
        rTag = new DataTag();
        dr = box.getSpace().makeVector();
        r0 = box.getSpace().makeVector();
        setInsertionDistances(new double[]{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1});
    }

    public void setInsertionDistances(double[] r) {
        rData = new DataDoubleArray(new int[]{r.length}, r);
        data = new DataFunction(new int[]{r.length});
        rDataInfo = new DataDoubleArray.DataInfoDoubleArray("r", Length.DIMENSION, new int[]{r.length});
        rDataInfo.addTag(rTag);
        dataInfo = new DataFunction.DataInfoFunction("y(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }

    /**
     * Sets the species, takes a prototype molecule, and gets handle to
     * appropriate species agent in box
     */
    public void setSpecies(ISpecies s) {
        species = s;
    }

    /**
     * Accessor for the species for which chemical potential is evaluated
     */
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * Number of Widom insertions attempted with each call to currentValue
     */
    public void setNInsert(int n) {
        nInsert = n;
    }

    /**
     * Accessor to number of Widom insertions attempted with each call to
     * currentValue
     */
    public int getNInsert() {
        return nInsert;
    }

    public IData getData() {
        double[] y = data.getData();
        IAtomList atoms = box.getLeafList();
        double[] r = rData.getData();
        for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            IAtom atom0 = atoms.get(random.nextInt(atoms.size()));
            r0.E(atom0.getPosition());
            dr.setRandomSphere(random);
            // move atom0 out of the way
            atom0.getPosition().PEa1Tv1(-1.5, dr);
            for (int j = 0; j < r.length; j++) {
                if (r[j] == 0) {
                    y[j]++;
                    continue;
                }
                int finalJ = j;
                IMolecule molecule = box.addNewMolecule(species, mol -> {
                    IAtom testAtom = mol.getChildList().get(0);
                    Vector testR = testAtom.getPosition();
                    testR.E(r0);
                    testR.PEa1Tv1(r[finalJ], dr);
                    Vector shift = box.getBoundary().centralImage(testR);
                    testR.PE(shift);
                    energyMeter.setTarget(testAtom);
                });
                double u = energyMeter.getDataAsScalar();
                box.removeMolecule(molecule);
                if (u == 0) y[j]++;
            }
            atom0.getPosition().E(r0);
        }
        for (int j = 0; j < r.length; j++) {
            y[j] /= nInsert;
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return rData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return rDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return rTag;
    }
}
