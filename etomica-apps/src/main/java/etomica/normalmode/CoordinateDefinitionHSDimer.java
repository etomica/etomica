/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.MoleculeChildAtomAction;
import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.IConformationOriented;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.models.nitrogen.AtomActionTransformed;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeArrayList;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * CoordinateDefinition implementation for HS Dimer. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position and 2 rotational displacements. 
 * 
 * @author Tai Boon Tan
 */
public class CoordinateDefinitionHSDimer extends CoordinateDefinitionMolecule
        implements Serializable {
	
    public CoordinateDefinitionHSDimer(Simulation sim, Box box, Primitive primitive, Basis basis, Space _space) {
    	
    	super(sim, box, primitive, 2, basis, _space);
    	
    	axis = space.makeVector();
    	
        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
        
        orientation = space.makeOrientation();
    }

    public void initializeCoordinates(int[] nCells) {
        IMoleculeList moleculeList = box.getMoleculeList();

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        Vector offset = lattice.getSpace().makeVector();
//        Vector[] primitiveVectors = primitive.vectors();
//        for (int i=0; i<primitiveVectors.length; i++) {
//            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
//        }
//        
//        offset.TE(-0.5);
        
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
        indexIterator.reset();
        Vector position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
		
        if (configuration != null){
        	configuration.initializeCoordinates(box);
        }

        //Scale the lattice offset
        Vector vectorOfMax = space.makeVector();
        Vector vectorOfMin = space.makeVector();
        Vector site = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);
        
        while (indexIterator.hasNext()) {
            site.E((Vector) lattice.site(indexIterator.next()));
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);
        
        indexIterator.reset();
                
        for (int iMolecule = 0; iMolecule<moleculeList.size(); iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IConformationOriented conformation = (IConformationOriented)(molecule.getType()).getConformation();
            
            if (configuration == null) {
            	// initialize the coordinate
                Vector v = space.makeVector();
                int thisOrientationIndex = orientationSelector.f(iMolecule);
                v.Ea1Tv1(cosTheta[thisOrientationIndex], axes[thisOrientationIndex][0]);
                v.PEa1Tv1(Math.sqrt(1-cosTheta[thisOrientationIndex]*cosTheta[thisOrientationIndex]), axes[thisOrientationIndex][1]);
                orientation.setDirection(v);

                conformation.initializePositions(molecule.getChildList(), orientation);
            }
                        
            int[] ii = indexIterator.next();
            // ii[0] and ii[1] = unit Cell number
            // ii[2] = molecule number in unit cell
        	//System.out.println(ii[0] +" " + ii[1] + " " + ii[2] + " " + ii[3] );
            
            position.E((Vector)lattice.site(ii));
            position.PE(offset);
            
            if (configuration == null) {
                atomActionTranslateTo.setDestination(position);
                atomActionTranslateTo.actionPerformed(molecule);
            }
            
            if (ii[space.D()] == 0) {
                // new cell
                iCell++;
                currentList = new MoleculeArrayList(basisSize);
                cells[iCell] = new BasisCell(currentList, lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }

            currentList.add(molecule);
        }

        moleculeSiteManager = new MoleculeAgentManager(sim.getSpeciesManager(), box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<>(new SiteSource(space), box);
    }
    
    public void setConfiguration(Configuration configuration){
        this.configuration = configuration;
    }

    public void setOrientations(Vector[][] newAxes, double[] newTheta, IntegerFunction newOrientationSelector){
        axes = newAxes;
        orientationSelector = newOrientationSelector;
        cosTheta = new double[newTheta.length];
        for (int i=0; i<axes.length; i++) {
            cosTheta[i] = Math.cos(newTheta[i]);
        }
    }
    
     public double[] calcU(IMoleculeList molecules) {
        
    	super.calcU(molecules);
    	
        int j = 3;
        
        for (int i = 0; i < molecules.size() ; i++){
        	IMolecule molecule = molecules.get(i);
        	int thisOrientationIndex = orientationSelector.f(i);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
	    	Vector leafPos1 = molecule.getChildList().get(1).getPosition();
	    	
	    	axis.Ev1Mv2(leafPos1, leafPos0);
	       	axis.normalize();
	       	
	       	double myCosTheta = axis.dot(axes[thisOrientationIndex][0]);
	       	if (myCosTheta < 0) {
	       	    throw new RuntimeException("oops");
	       	}
            u[j] = myCosTheta - cosTheta[thisOrientationIndex];

            double sinthetaOver2 = 0.5*Math.sqrt(axis.Mv1Squared(axes[thisOrientationIndex][0]));
	       	double sintheta = 2*Math.sqrt(1-sinthetaOver2*sinthetaOver2)*sinthetaOver2;
	    	
	       	double x = axis.dot(axes[thisOrientationIndex][1]);  
	    	double y = axis.dot(axes[thisOrientationIndex][2]);
	    	double phi = Math.atan2(y, x);
	    	
	    	if (sintheta == 0){
	    	    // not good!
	    	    u[j+1] = 0;
	    	}
	    	else {
	    	    u[j+1] = phi/sintheta;
	    	}
	    	j += coordinateDim/molecules.size();
        }
        return u;
     }

     /**
      * return the initial Orientation of the molecule
      */
    public IOrientation getMoleculeOrientation(IMolecule molecule) {
        Vector v = space.makeVector();
        int thisOrientationIndex = orientationSelector.f(molecule.getIndex());
        v.Ea1Tv1(cosTheta[thisOrientationIndex], axes[thisOrientationIndex][0]);
        v.PEa1Tv1(Math.sqrt(1-cosTheta[thisOrientationIndex]*cosTheta[thisOrientationIndex]), axes[thisOrientationIndex][1]);
        orientation.setDirection(v);
        return orientation;
    }
    
    public void setToU(IMoleculeList molecules, double[] newU) {
    	
    	int j=3;
	        
        for (int i = 0; i < molecules.size() ; i++){
        	
        	IMolecule molecule = molecules.get(i);
            int thisOrientationIndex = orientationSelector.f(i);
            Vector[] myAxes = axes[thisOrientationIndex];

            double myCosTheta = cosTheta[thisOrientationIndex] + newU[j];
            double mySinTheta = Math.sqrt(1-myCosTheta*myCosTheta);
            double myPhi = newU[j+1] / mySinTheta;
            double mySinPhi = Math.sin(myPhi);
            double myCosPhi = Math.cos(myPhi);
            
            axis.Ea1Tv1(myCosTheta, myAxes[0]);
            axis.PEa1Tv1(mySinTheta*myCosPhi, myAxes[1]);
            axis.PEa1Tv1(mySinTheta*mySinPhi, myAxes[2]);
            orientation.setDirection(axis);

            IConformationOriented conformationOriented = (IConformationOriented) molecule.getType().getConformation();
            conformationOriented.initializePositions(molecule.getChildList(), orientation);
	    	
	    	j += coordinateDim/molecules.size();
	    	
        }
        super.setToU(molecules, newU);
    }

    
    public void setToUMoleculei(int moleculei, double[] newU) {
    	
    	if(newU.length != 5){
    		throw new RuntimeException("<CoordinateDefinitionNitrogen> setToUMoleculei method, newU[] length should be 5!");
    	}
    	
		IMolecule molecule = box.getMoleculeList().get(moleculei);
        int thisOrientationIndex = orientationSelector.f(moleculei);
        Vector[] myAxes = axes[thisOrientationIndex];

        double myCosTheta = cosTheta[thisOrientationIndex] + newU[3];
        double mySinTheta = Math.sqrt(1-myCosTheta*myCosTheta);
        double myPhi = newU[4] / mySinTheta;
        double mySinPhi = Math.sin(myPhi);
        double myCosPhi = Math.cos(myPhi);
        
        axis.Ea1Tv1(myCosTheta, myAxes[0]);
        axis.PEa1Tv1(mySinTheta*myCosPhi, myAxes[1]);
        axis.PEa1Tv1(mySinTheta*mySinPhi, myAxes[2]);
        orientation.setDirection(axis);

        IConformationOriented conformationOriented = (IConformationOriented) molecule.getType().getConformation();
        conformationOriented.initializePositions(molecule.getChildList(), orientation);
		        
        Vector site = getLatticePosition(molecule);
        for (int k = 0; k < site.getD(); k++) {
        	work1.setX(k, site.getX(k) + newU[k]);
        }
		            
        atomActionTranslateTo.setDestination(work1);
        atomActionTranslateTo.actionPerformed(molecule);
    }
    
    private static final long serialVersionUID = 1L;

    protected Vector[][] axes;
    protected final IOrientation orientation;
    protected IntegerFunction orientationSelector;
    protected final Vector axis;
    protected Configuration configuration;
    protected final MoleculeChildAtomAction atomGroupAction;
    protected double[] cosTheta;

    public interface IntegerFunction {
        public int f(int i);
    }
}
