/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

/**
 * Creates a configuration using a CubicLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLatticeTube extends ConfigurationLattice {

    public ConfigurationLatticeTube(BravaisLatticeCrystal lattice,
    		               double length, Space _space) {
        this(lattice, length, new IndexIteratorRectangular(lattice.D()), _space);//need a default iterator
    }
	/**
	 * Constructor for ConfigurationLatticeTube.
	 * @param space
	 */
	public ConfigurationLatticeTube(BravaisLatticeCrystal lattice,
			        double length, IndexIteratorSizable indexIterator,
			        Space _space) {
	    super(lattice, _space);
        this.indexIterator = indexIterator;
        this.length = length;
        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
	}
	
	public void setSpeciesSpheres(SpeciesSpheresMono[] speciesSpheres) {
	    this.speciesSpheres = speciesSpheres;
	}
	
	public void setSpeciesTube(SpeciesSpheres speciesTube) {
	    this.speciesTube = speciesTube;
	}
	
    public void initializeCoordinates(Box box) {
        IMoleculeList[] spheresLists = new IMoleculeList[]{box.getMoleculeList(speciesSpheres[0]), box.getMoleculeList(speciesSpheres[1])};
        
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int)Math.ceil((double)spheresLists[0].size()/(double)basisSize);
        
        //determine scaled shape of simulation volume
        Vector shape = space.makeVector();
        shape.E(box.getBoundary().getBoxSize());
        shape.setX(2,shape.getX(2)*length);
        Vector latticeConstantV = Vector.of(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = basisSize;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }
    
        // determine lattice constant
        Vector latticeScaling = space.makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(shape);
            latticeScaling.DE(Vector.of(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        Vector offset = space.makeVector();
        offset.E(box.getBoundary().getBoxSize());
        Vector vectorOfMax = space.makeVector();
        Vector vectorOfMin = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        // XXX this can do strange things. it's probably not needed for 
        // periodic boundaries, but gets the atoms off the boundaries for 
        // non-periodic boundaries
        indexIterator.reset();
        while (indexIterator.hasNext()) {
            Vector site = (Vector) lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);
        offset.setX(2, offset.getX(2) - 0.5*box.getBoundary().getBoxSize().getX(2)*(1-length));

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules  
        indexIterator.reset();
        
        // first species (mono spheres)
        int nSpheres = spheresLists[0].size();
        for (int i=0; i<nSpheres; i++) {
            IMolecule a = spheresLists[0].get(i);
            
            int[] ii = indexIterator.next();
            Vector site = (Vector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        double z = offset.getX(2);
        offset.setX(2,z+box.getBoundary().getBoxSize().getX(2)*(1-length));
        myLat = new MyLattice(lattice, latticeScaling, offset);
        indexIterator.reset();
        
        nSpheres = spheresLists[1].size();
        // second species (mono spheres)
        for (int i=0; i<nSpheres; i++) {
            IMolecule a = spheresLists[1].get(i);
            
            int[] ii = indexIterator.next();
            Vector site = (Vector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        //loop for multiple tubes.
        IMoleculeList tubeList = box.getMoleculeList(speciesTube);
        int nTubes = tubeList.size();
        atomActionTranslateTo.setAtomPositionDefinition(new MoleculePositionGeometricCenter(space));
        // put them all at 0.  oops
        atomActionTranslateTo.setDestination(space.makeVector());
        for (int i=0; i<nTubes; i++) {
            IMolecule a = tubeList.get(i);
            a.getType().initializeConformation(a);
            atomActionTranslateTo.actionPerformed(a);
        }
        
    }
    
    private static final long serialVersionUID = 1L;
    private final IndexIteratorSizable indexIterator;
    private final MoleculeActionTranslateTo atomActionTranslateTo;
    protected SpeciesSpheresMono[] speciesSpheres;
    protected SpeciesSpheres speciesTube;
    private final double length;
}
