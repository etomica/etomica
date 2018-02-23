/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.box.Box;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;

/**
 * Sets the configuration to the zincblende structure, which consists
 * of two fcc lattices, with one shifted in each direction by one-quarter
 * of the lattice constant.
 */
public class ConfigurationZincblende extends ConfigurationLattice {
    
    private static final long serialVersionUID = 2L;
    private MoleculeChildAtomAction translator0, translator1;
    protected ISpecies[] species;
    
    public ConfigurationZincblende(double latticeConstant, Space space) {
        super(new LatticeCubicFcc(space, latticeConstant), space);
        species = new ISpecies[2];
    }
    
    public void setSpecies1(ISpecies species1) {
        species[0] = species1;
    }
    
    public void setSpecies2(ISpecies species2) {
        species[1] = species2;
    }
    
    public ISpecies getSpecies1() {
        return species[0];
    }
    
    public ISpecies getSpecies2() {
        return species[1];
    }
    
    /**
     * Initializes positions of atoms to the zincblende structure.  The given
     * array should hold exactly two AtomLists, each with the same number of atoms.
     */
    public void initializeCoordinates(Box box) {
        translator0 = new MoleculeChildAtomAction(new AtomActionTranslateBy(space));
        translator1 = new MoleculeChildAtomAction(new AtomActionTranslateBy(space));
        IMoleculeList[] lists = new IMoleculeList[]{box.getMoleculeList(species[0]), box.getMoleculeList(species[1])};
        if(lists == null || lists.length != 2) {//need an exception for this
            throw new IllegalArgumentException("inappropriate argument to ConfigurationZincBlende");
        }
        if(lists[0].size() != lists[1].size()) {
            System.err.println("Warning: different numbers of molecules for two species in ConfigurationZincBlende");
        }
        
        int nCells = (int) Math.ceil(lists[0].size() / 4.0);

        // determine scaled shape of simulation volume
        Vector shape = space.makeVector();
        shape.E(box.getBoundary().getBoxSize());
        Vector latticeConstantV = space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = 4;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        //shift lattice in all three directions by one-quarter the lattice constant
        Vector3D shift = new Vector3D();
        shift.Ea1Tv1(-0.5,box.getBoundary().getBoxSize());
        shift.PE(0.125*((LatticeCubicFcc)lattice).getLatticeConstant());
        ((AtomActionTranslateBy)translator0.getAtomAction()).setTranslationVector(shift);

        shift.PE(0.25*((LatticeCubicFcc)lattice).getLatticeConstant());
        ((AtomActionTranslateBy)translator1.getAtomAction()).setTranslationVector(shift);

        // Place molecules
        indexIterator.reset();
        int i = 0;
        while (indexIterator.hasNext()) {
            int[] ii = indexIterator.next();
            Vector site = (Vector) lattice.site(ii);
            atomActionTranslateTo.setDestination(site);

            IMolecule a0 = lists[0].get(i);
            IMolecule a1 = lists[1].get(i);
            atomActionTranslateTo.actionPerformed(a0);
            atomActionTranslateTo.actionPerformed(a1);

            translator0.actionPerformed(a0);
            translator1.actionPerformed(a1);
        }
    }        
}
