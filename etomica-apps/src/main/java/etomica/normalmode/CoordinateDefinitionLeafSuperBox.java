/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.ISpecies;
import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.MoleculeListWrapper;
import etomica.molecule.iterator.MoleculeIteratorAllMolecules;
import etomica.space.Space;
import etomica.space.Vector;



public class CoordinateDefinitionLeafSuperBox extends CoordinateDefinitionLeaf {


	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public CoordinateDefinitionLeafSuperBox(Box box, Primitive primitive,
                                            Basis basis, Space space) {
		super(box, primitive, basis, space);
	}
	
    public void initializeCoordinates(int[] nCells) {
        MoleculeIteratorAllMolecules atomIterator = new MoleculeIteratorAllMolecules(box);
        IMoleculeList moleculeList = box.getMoleculeList();
        if (moleculeList.getMoleculeCount() == 0) {
            throw new RuntimeException("There are no atoms yet!");
        }

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        Vector offset = lattice.getSpace().makeVector();
        Vector[] primitiveVectors = primitive.vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
        }
        offset.TE(-0.5);
        
        IndexIteratorRectangular indexIterator = new IndexIteratorRectangular(space.D()+1);
        int[] iteratorDimensions = new int[space.D()+1];
        
        
        System.arraycopy(nCells, 0, iteratorDimensions, 0, nCells.length);
        iteratorDimensions[nCells.length] = basisSize;
        indexIterator.setSize(iteratorDimensions);

        int totalCells = 1;
        for (int i=0; i<nCells.length; i++) {
            totalCells *= nCells[i];
        }
        
        cells = new BasisCell[totalCells];
        int iCell = -1;
        // Place molecules
        atomIterator.reset();
        indexIterator.reset();
        Vector position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
        
        int counterSpeciesA =0;
        int counterSpeciesB =0;
        
        for (int[] ii = indexIterator.next(); ii != null; ii = indexIterator.next()){
        	
        	boolean inCenterBox = true;
        	
        	if(is864){
        		for (int i=0; i<space.D(); i++){
        			if (ii[i] < nCells[i]/3 || ii[i] > (nCells[i]*2/3-1)){ //outside the box
        				inCenterBox = false;
        			
        			}
        		}
        	} else {
        		for (int i=0; i<space.D(); i++){
        			if (ii[i] < nCells[i]/4 || ii[i] > (nCells[i]*3/4-1)){ //outside the box
        				inCenterBox = false;
        			
        			}
        		}
        	}
        	
        	
        	IMolecule molecule;
        	if (inCenterBox) {
        		molecule = box.getMoleculeList(speciesA).getMolecule(counterSpeciesA);
        		counterSpeciesA ++;
        		
        	} else {
        		molecule = box.getMoleculeList(speciesB).getMolecule(counterSpeciesB);
        		counterSpeciesB ++;
        		
        	}
        	      	       	
        	// initialize coordinates of child atoms
        	molecule.getType().initializeConformation(molecule);

            position.E((Vector)lattice.site(ii));
            position.PE(offset);
            
            atomActionTranslateTo.setDestination(position);
            atomActionTranslateTo.actionPerformed(molecule);

            if (ii[space.D()] == 0) {
                if (iCell > -1) {
                    initNominalU(cells[iCell].molecules);
                }
                // new cell
                iCell++;
                currentList = new MoleculeArrayList(basisSize);
                cells[iCell] = new BasisCell(new MoleculeListWrapper(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(molecule);

        }
                
        initNominalU(cells[totalCells-1].molecules);
        
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box, Vector.class);
    }
    
    public void setSpecies(ISpecies speciesA, ISpecies speciesB){
    	this.speciesA = speciesA;
    	this.speciesB = speciesB;
    }
    
    public void setIs864(){
    	is864 = true;
    }
    
    public void setIs256(){
    	is864 = false;
    }
    
    public boolean is864(){
    	return is864;
    }
    
    protected ISpecies speciesA, speciesB; 
    protected boolean is864 = true;
    
}
