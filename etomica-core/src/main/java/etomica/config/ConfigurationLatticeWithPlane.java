/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.box.Box;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.SpaceLattice;
import etomica.math.geometry.Plane;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Constructs configuration that has the molecules placed on the sites of a
 * lattice. Molecules are placed sequentially on the sites of the lattice in an
 * order specified by an index iterator set at construction.
 * <p>
 * Iteration over the indexes yields integer arrays, and with each iteration the
 * array is passed to the site method of the lattice which returns the position
 * for placement of the next molecule. Each array index is iterated to a maximum
 * value determined by the number of molecules to be placed, the dimensions of
 * the box in which they are placed, and the lattice constants of the lattice.
 * <p>
 * An instance of this class may be configured to place atoms such that they
 * uniformly fill the volume of the box. It will attempt this by scaling the
 * lattice constants of the configuration in an appropriate way. Success in
 * getting a good spatial distribution may vary.
 * <p>
 * An instance can also be configured to remember the indices used to get
 * lattice position for each molecule that is placed. This can be useful if it
 * is desired to associate each molecule with a lattice site.
 */
public class ConfigurationLatticeWithPlane extends ConfigurationLattice {

    private Plane plane;
    private ArrayList<ISpecies> species;
    private HashMap<ISpecies,Float> allocation;
    
    private final int LEFT  = 0;
    private final int RIGHT = 1;

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationLatticeWithPlane(SpaceLattice lattice, Plane plane, Space space) {
        this(lattice, plane, new IndexIteratorRectangular(lattice.D()), space);
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    private ConfigurationLatticeWithPlane(SpaceLattice lattice,
            Plane plane, IndexIteratorSizable indexIterator, Space space) {
    	super(lattice, indexIterator, space);

        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.plane = plane;
        species = new ArrayList<ISpecies>();
        allocation = new HashMap<ISpecies,Float>();
    }

    /**
     * 
     * @param newSpecies : species that will need to be configured
     */
    public void addSpecies(ISpecies newSpecies) {
    	if(species.indexOf(newSpecies) == -1) {
    	    species.add(newSpecies);
    	    allocation.put(newSpecies, 0.5F);
    	}
    }

    /**
     * 
     * @param remSpecies : species that should be removed from configuration
     */
    public void removeSpecies(ISpecies remSpecies) {
    	if(species.indexOf(remSpecies) != -1) {
    		species.remove(remSpecies);
    	    allocation.remove(remSpecies);
    	}
    }

    /**
     * Sets the initial allocation of the species on the left side
     * of the plane.  The remaining molecules in the species will
     * be initially allocated to the right side of the plane.
     * 
     * @param sp : Species to set allocation for
     * @param pct : Percentage of species to initially allocate
     * on the left side of the plane (0.0 <= pct <= 1.0)
     */
    public void setSpeciesAllocation(ISpecies sp, float pct) {
    	if(species.indexOf(sp) != -1 && pct <= 1.0) {
    		allocation.remove(sp);
    		allocation.put(sp, pct);
    	}
    }

    /**
     * Sets the initial allocation of the species on the left side
     * of the plane.  The remaining molecules in the species will
     * be initially allocated to the right side of the plane.
     * 
     * @param sp : Species to set allocation for
     */
    public float getSpeciesAllocation(ISpecies sp) {
    	if(allocation.containsKey(sp)) {
    		return allocation.get(sp);
     	}
		return 0;
    }

    /**
     * Places the molecules in the given box on the positions of the
     * lattice.  The molecules must be of a species added with addSpecies()
     * or they will not be displayed.
     */
// Method is assuming plane is in the middle.  Calculations making this
// assumption are noted.
    public void initializeCoordinates(Box box) {

        int numSpecies = species.size();
        int[] speciesCount = new int[numSpecies];

        int sumOfMolecules = 0;
        for (int i = 0; i < numSpecies; i++) {
            speciesCount[i] = box.getNMolecules(species.get(i));
            sumOfMolecules = sumOfMolecules + speciesCount[i];
        }

        if (sumOfMolecules == 0) {
            return;
        }

        // Allocate species to left or right side of plane
        int [][] molecules = new int[RIGHT+1][numSpecies];
        int [] maxMolecules = { 0, 0 };

        for (int i = 0; i < numSpecies; i++) {
        	ISpecies sp = species.get(i);
            molecules[LEFT][i] =  (int)(box.getNMolecules(sp) * getSpeciesAllocation(sp));
            molecules[RIGHT][i] = box.getNMolecules(sp) - molecules[LEFT][i];
            maxMolecules[LEFT] += molecules[LEFT][i];
            maxMolecules[RIGHT] += molecules[RIGHT][i];
        }

        // determine scaled shape of simulation volume
        Vector halfShape = space.makeVector();
        halfShape.E(box.getBoundary().getBoxSize());

	    int planeDimIdx = 0;
	    if(plane.getA() > plane.epsilon) {
	        planeDimIdx = 0;
	    }
	    else if(plane.getB() > plane.epsilon) {
	    	planeDimIdx = 1;
	    }
	    else if(space.D() == 3 &&
	    		plane.getC() > plane.epsilon) {
	    	planeDimIdx = 2;
	    }

	    Vector entireShape = space.makeVector();
	    entireShape.E(halfShape);

	    //  NOTE, JUST DIVIDING BY 2 ASSUMES PLANE DOWN CENTER
        // Need to adjust shape on either side of plane in plane dimension.
	    halfShape.setX(planeDimIdx, halfShape.getX(planeDimIdx) / 2);


        Vector latticeConstantV = Vector.of(lattice.getLatticeConstants());
        halfShape.DE(latticeConstantV);

        int[][] latticeDimensions;
        latticeDimensions = new int[RIGHT+1][];

        for(int side = LEFT; side <= RIGHT; side++) {
		    // determine number of cells in each direction
		    latticeDimensions[side] = calculateLatticeDimensions(maxMolecules[side], halfShape);
        }

        int[] molIdx = new int[numSpecies];
        for(int side = LEFT; side <= RIGHT; side++) {

	        // determine lattice constant
	        Vector latticeScaling = space.makeVector();
	        if (rescalingToFitVolume) {
                latticeScaling.E(halfShape);
	            latticeScaling.DE(Vector.of(latticeDimensions[side]));
	        } else {
	            latticeScaling.E(1.0);
	        }

            indexIterator.setSize(latticeDimensions[side]);
            indexIterator.reset();

	        // determine amount to shift lattice so it is centered in volume
	        Vector offset = space.makeVector();
	        offset.E(box.getBoundary().getBoxSize());

	        Vector temp3 = space.makeVector();
            temp3.E(entireShape);
            temp3.TE(-0.5);
            temp3.setX(planeDimIdx, halfShape.getX(planeDimIdx) * (side-1));
            offset.E(temp3);

	        Vector temp2 = space.makeVector();
	        temp2.E(latticeScaling);
	        temp2.TE(0.5);
	        offset.PE(temp2);

	        myLat = new MyLattice(lattice, latticeScaling, offset);
	
	        indexIterator.reset();

	        // Place molecules
	        for (int i = 0; i < numSpecies; i++) {
                IMoleculeList imolecules = box.getMoleculeList(species.get(i));
	        	for (int y = 0; y < molecules[side][i]; y++) {
	                IMolecule a = imolecules.get(molIdx[i]);
                    molIdx[i]++;

			        int[] idx = indexIterator.next();

			        atomActionTranslateTo.setDestination(myLat.site(idx));
			        atomActionTranslateTo.actionPerformed(a);
	
		        }
	        }
        }

    }
}
