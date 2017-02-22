/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.Serializable;

import etomica.action.BoxInflate;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.space.ISpace;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz
 */
public class CoordinateDefinitionMolecule extends CoordinateDefinition
        implements Serializable {

    public CoordinateDefinitionMolecule(ISimulation sim, IBox box, Primitive primitive, int orientationDim, ISpace space) {
        this(sim, box, primitive, orientationDim, new BasisMonatomic(space), space);
    }
    
    public CoordinateDefinitionMolecule(ISimulation sim, IBox box, Primitive primitive, int orientationDim, Basis basis, ISpace space) {
        super(box, (space.D() + orientationDim)*basis.getScaledCoordinates().length, primitive, basis, space);
        this.sim = sim;
        work1 = space.makeVector();
        inflate = new BoxInflate(space);
        inflate.setBox(box);
       
        u = new double[coordinateDim];
        setPositionDefinition(new AtomPositionGeometricCenter(space));
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
    
    	int j = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IVector pos = positionDefinition.position(molecule);
            IVectorMutable site = getLatticePosition(molecule);
            
            work1.Ev1Mv2(pos, site);
               
            for (int k = 0; k < pos.getD(); k++) {
                u[j+k] = work1.getX(k);
            }
               j += coordinateDim/molecules.getMoleculeCount();

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
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IVectorMutable site = getLatticePosition(molecule);
            for (int k = 0; k < site.getD(); k++) {
                work1.setX(k, site.getX(k) + newU[j+k]);
            }
            
            atomActionTranslateTo.setDestination(work1);
            atomActionTranslateTo.actionPerformed(molecule);
            
            j += coordinateDim/molecules.getMoleculeCount();

        }
    }
    
    public IVectorMutable getLatticePosition(IMolecule molecule) {
        return (IVectorMutable)moleculeSiteManager.getAgent(molecule);
    }
    
    public void setPositionDefinition(IAtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
        atomActionTranslateTo.setAtomPositionDefinition(positionDefinition);
    }

    public IAtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }
    
    public void setInitVolume(IVector initV){
    	this.initVolume = initV;
    }
    
    private static final long serialVersionUID = 1L;
    protected final ISimulation sim;
    protected MoleculeAgentManager moleculeSiteManager;
    protected final IVectorMutable work1;
    protected final double[] u;
    protected IAtomPositionDefinition positionDefinition;
    protected double rScale;
    protected IVector initVolume;
    protected final BoxInflate inflate;

    protected static class MoleculeSiteSource implements MoleculeAgentSource, Serializable {
        
        public MoleculeSiteSource(ISpace space, IAtomPositionDefinition positionDefinition) {
            this.space = space;
            this.positionDefinition = positionDefinition;
        }
        public Class getMoleculeAgentClass() {
            return IVectorMutable.class;
        }
        public Object makeAgent(IMolecule molecule) {
            IVectorMutable vector = space.makeVector();
            vector.E(positionDefinition.position(molecule));
            return vector;
        }
        public void releaseAgent(Object agent, IMolecule molecule) {
            //nothing to do
        }

        private final ISpace space;
        protected final IAtomPositionDefinition positionDefinition;
        private static final long serialVersionUID = 1L;
    }
}
