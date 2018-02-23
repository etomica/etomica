/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.box.Box;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.*;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * with an extra degree of freedom for volume fluctuation
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class CoordinateDefinitionMoleculeVolumeFluctuation extends CoordinateDefinition
        implements Serializable {

    public CoordinateDefinitionMoleculeVolumeFluctuation(Simulation sim, Box box, Primitive primitive, int orientationDim, Space space) {
        this(sim, box, primitive, orientationDim, new BasisMonatomic(space), space);
    }
    
    public CoordinateDefinitionMoleculeVolumeFluctuation(Simulation sim, Box box, Primitive primitive, int orientationDim, Basis basis, Space space) {
        super(box, ((space.D() + orientationDim)*basis.getScaledCoordinates().length+1), primitive, basis, space);
        this.sim = sim;
        work1 = space.makeVector();
        inflate = new BoxInflate(space);
        inflate.setBox(box);
        
        u = new double[coordinateDim];
        setPositionDefinition(new MoleculePositionGeometricCenter(space));
        rScale = 1.0;
    }
    
    public void initializeCoordinates(int[] nCells) {
        super.initializeCoordinates(nCells);
        moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
    }

    public double[] calcU(IMoleculeList molecules) {
        // calculates components of U related to the the center of mass of the
        // molecules
        // subclass is responsible for setting orientation or intramolecular
        // degrees of freedom
    	
    	/*
    	 * assigning the degrees of freedom for volume fluctuation
    	 * in the last 3 components in u[] array
    	 * 
    	 *  x-direction fluctuation : u[coordinateDim-3]
    	 *  y-direction fluctuation : u[coordinateDim-2]
    	 *  z-direction fluctuation : u[coordinateDim-1]
    	 *  
    	 */
    	
    	double currentDimi = box.getBoundary().getBoxSize().getX(0);
    	rScale = currentDimi/initVolume.getX(0);
    	u[coordinateDim-1] = currentDimi - initVolume.getX(0);  
    		
    	
        
    	int j = 0;
        for (int i = 0; i<molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            Vector pos = positionDefinition.position(molecule);
            Vector site = getLatticePosition(molecule);
            
            work1.E(new double[]{(1/rScale)*pos.getX(0),
            					 (1/rScale)*pos.getX(1),
            					 (1/rScale)*pos.getX(2)});
            
            work1.ME(site);
            //work1.Ev1Mv2(pos, site);
               
            for (int k = 0; k < pos.getD(); k++) {
                u[j+k] = work1.getX(k);
            }
               j += (coordinateDim-1)/molecules.size();

        }
        
        return u;
    }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(IMoleculeList molecules) {
    }

    public void setToU(IMoleculeList molecules, double[] newU) {
        // sets the center of mass of the molecules to that specified by newU
        // subclass is responsible for setting orientation or intramolecular
        // degrees of freedom
    	

    	
        int j = 0;
        for (int i = 0; i<molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            Vector site = getLatticePosition(molecule);
            for (int k = 0; k < site.getD(); k++) {
                work1.setX(k, site.getX(k) + newU[j+k]);
            }
            
            atomActionTranslateTo.setDestination(work1);
            atomActionTranslateTo.actionPerformed(molecule);
            
            j += coordinateDim/molecules.size();

        }
    	
        rScale = (initVolume.getX(0) + newU[coordinateDim-1])/initVolume.getX(0); //rescale the fluctation to assigned fluctuation
    	
       	inflate.setScale(rScale);
    	inflate.actionPerformed();
    	    	   	
    }
    
    public Vector getLatticePosition(IMolecule molecule) {
        return (Vector)moleculeSiteManager.getAgent(molecule);
    }
    
    public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
        atomActionTranslateTo.setAtomPositionDefinition(positionDefinition);
    }

    public IMoleculePositionDefinition getPositionDefinition() {
        return positionDefinition;
    }
    
    public void setInitVolume(Vector initV){
    	this.initVolume = initV;
    }
    
    private static final long serialVersionUID = 1L;
    protected final Simulation sim;
    protected MoleculeAgentManager moleculeSiteManager;
    protected final Vector work1;
    protected final double[] u;
    protected IMoleculePositionDefinition positionDefinition;
    protected double rScale;
    protected Vector initVolume;
    protected final BoxInflate inflate;

    protected static class MoleculeSiteSource implements MoleculeAgentSource, Serializable {
        
        public MoleculeSiteSource(Space space, IMoleculePositionDefinition positionDefinition) {
            this.space = space;
            this.positionDefinition = positionDefinition;
        }

        public Object makeAgent(IMolecule molecule) {
            Vector vector = space.makeVector();
            vector.E(positionDefinition.position(molecule));
            return vector;
        }
        public void releaseAgent(Object agent, IMolecule molecule) {
            //nothing to do
        }

        private final Space space;
        protected final IMoleculePositionDefinition positionDefinition;
        private static final long serialVersionUID = 1L;
    }
}
