package etomica.paracetamol;

import java.io.Serializable;

import etomica.action.AtomGroupAction;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomsetArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Conformation;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinitionMolecule;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Tensor;
import etomica.space3d.IVector3D;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinitionMolecule
        implements Serializable {

    public CoordinateDefinitionParacetamol(Phase phase, Primitive primitive, Basis basis) {
    	super(phase, primitive, 3, basis);
       
       	axes = new IVector3D [3];
        axes [0] = (IVector3D)phase.getSpace().makeVector();
        axes [1] = (IVector3D)phase.getSpace().makeVector();
        axes [2] = (IVector3D)phase.getSpace().makeVector();
        com = (IVector3D)phase.getSpace().makeVector();
        temp = (IVector3D)phase.getSpace().makeVector();
        proj = (IVector3D)phase.getSpace().makeVector();
        proja = (IVector3D)phase.getSpace().makeVector();
        projb = (IVector3D)phase.getSpace().makeVector();
        axisNorm = (IVector3D)phase.getSpace().makeVector();
        axisNormPrime = (IVector3D)phase.getSpace().makeVector();
        axis0 = (IVector3D)phase.getSpace().makeVector();
        axis0Prime = (IVector3D)phase.getSpace().makeVector();
        b = (IVector3D)phase.getSpace().makeVector();
        bprime = (IVector3D)phase.getSpace().makeVector();
        c = (IVector3D)phase.getSpace().makeVector();
        vector1 = (IVector3D)phase.getSpace().makeVector();
        vector2 = (IVector3D)phase.getSpace().makeVector();
        Normal = (IVector3D)phase.getSpace().makeVector();
        atomGroupAction = new AtomGroupAction(new AtomActionTransformed(lattice.getSpace()));
    }

    public void initializeCoordinates(int[] nCells) {
        AtomIteratorAllMolecules atomIterator = new AtomIteratorAllMolecules(phase);

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        IVector offset = lattice.getSpace().makeVector();
        IVector[] primitiveVectors = primitive.vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
        }
        offset.TE(-0.5);
        
        IndexIteratorRectangular indexIterator = new IndexIteratorRectangular(phase.getSpace().D()+1);
        int[] iteratorDimensions = new int[phase.getSpace().D()+1];
        
        
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
        IVector position = lattice.getSpace().makeVector();
        AtomArrayList currentList = null;
        Tensor t = lattice.getSpace().makeTensor();
    	
        for (IAtom a = atomIterator.nextAtom(); a != null;
             a = atomIterator.nextAtom()) {
            if (a instanceof IAtomGroup) {
                // initialize coordinates of child atoms
                Conformation config = a.getType().creator().getConformation();
                config.initializePositions(((IAtomGroup)a).getChildList());
            }
            
            
            int[] ii = indexIterator.next();
            
            for (int i=0; i<3; i++){
            	t.setComponent(i, i, basisOrientation[ii[3]][i]);
            	
          	  ((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(t);
              atomGroupAction.actionPerformed(a);
            }
            
            position.E((IVector)lattice.site(ii));
            position.PE(offset);
            
            atomActionTranslateTo.setDestination(position);
            atomActionTranslateTo.actionPerformed(a);

            if (ii[phase.getSpace().D()] == 0) {
                if (iCell > -1) {
                    initNominalU(cells[iCell].molecules);
                }
                // new cell
                iCell++;
                currentList = new AtomArrayList(basisSize);
                cells[iCell] = new BasisCell(new AtomsetArrayList(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(a);
        }
        
        initNominalU(cells[totalCells-1].molecules);
        
        siteManager = new AtomAgentManager(new SiteSource(phase.getSpace()), phase);
    }
    
    public double [][] setBasisOrthorhombic(){
    	basisOrientation = new double [8][3];
    	
    	basisOrientation[0][0] =  1;
    	basisOrientation[0][1] =  1;
    	basisOrientation[0][2] =  1;
    	
    	basisOrientation[1][0] = -1;
    	basisOrientation[1][1] = -1;
    	basisOrientation[1][2] =  1;
    	
    	basisOrientation[2][0] = -1;
    	basisOrientation[2][1] =  1;
    	basisOrientation[2][2] = -1;
    	
    	basisOrientation[3][0] =  1;
    	basisOrientation[3][1] = -1;
    	basisOrientation[3][2] = -1;
    	
    	basisOrientation[4][0] = -1;
    	basisOrientation[4][1] = -1;
    	basisOrientation[4][2] = -1;
    	
    	basisOrientation[5][0] =  1;
    	basisOrientation[5][1] =  1;
    	basisOrientation[5][2] = -1;
    	
    	basisOrientation[6][0] =  1;
    	basisOrientation[6][1] = -1;
    	basisOrientation[6][2] =  1;
    	
    	basisOrientation[7][0] = -1;
    	basisOrientation[7][1] =  1;
    	basisOrientation[7][2] =  1;
    	
    	return basisOrientation;
    }
    
    public double [][] setBasisMonoclinic(){
    	basisOrientation = new double [4][3];
    	
    	basisOrientation[0][0] =  1;
    	basisOrientation[0][1] =  1;
    	basisOrientation[0][2] =  1;
    	
    	basisOrientation[1][0] = -1;
    	basisOrientation[1][1] =  1;
    	basisOrientation[1][2] = -1;
    	
    	basisOrientation[2][0] = -1;
    	basisOrientation[2][1] = -1;
    	basisOrientation[2][2] = -1;
    	
    	basisOrientation[3][0] =  1;
    	basisOrientation[3][1] = -1;
    	basisOrientation[3][2] =  1;
    	
    	return basisOrientation;
    }
    
    /*
     * 
     */
    
    public double[] calcU(AtomSet molecules) {
        super.calcU(molecules);
        
        int j = 3;
        
        for (int i=0; i < molecules.getAtomCount() ; i++){
        	
        	IAtomGroup molecule = (IAtomGroup)molecules.getAtom(i);
        	IVector3D [] siteOrientation = (IVector3D [])orientationManager.getAgent(molecule);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().get(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().get(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().get(10)).getPosition();
	    	
	    	/*
	    	 * Determine axis 1 by using Vector Projection
	    	 */
	    	axis0Prime.Ev1Mv2(leafPos5, leafPos0);
	    	
	    	axisNormPrime.Ev1Mv2(leafPos10, leafPos0);
	    	proj.E(axis0Prime);
	    	double dotProd = proj.dot(axisNormPrime);
	    	proj.Ea1Tv1(dotProd/axis0Prime.squared(), axis0Prime);
	    	axisNormPrime.ME(proj);
	    	
	      	axis0Prime.TE(1 /Math.sqrt(axis0Prime.squared()));
	      	axisNormPrime.TE(1 /Math.sqrt(axisNormPrime.squared()));
	    
	    	/*
	    	 * Determining the angle between axes[0] and siteOrientation[1,2]
	    	 * we have axes [0] and siteOrientation [0,1,2]
	    	 */ 
	    	
	    	u[j] = axis0Prime.dot(siteOrientation[1]);
	    	u[j+1] = axis0Prime.dot(siteOrientation[2]);
	      	
	    	/*
	    	 * Getting rotation angle at x-axis
	    	 * 
	    	 * axisNorm is perpendicular to siteOrientation[0]
	    	 */
        	double cosAlpha = axis0Prime.dot(siteOrientation[0]);
	    	if (cosAlpha > 0.99999999){
	    		axisNorm.E(axisNormPrime);
	    		
	    	} else {
	    		
	    		bprime.Ea1Tv1( 1.0 /cosAlpha, siteOrientation[0]);
	    		bprime.ME(axis0Prime);
	    		bprime.TE(1.0 /Math.sqrt(bprime.squared()));
	    		
	    		b.Ea1Tv1(cosAlpha, bprime);
	    		b.PEa1Tv1(-Math.sqrt(1-cosAlpha*cosAlpha), axis0Prime);
	    		
	    		double bComponent = axisNormPrime.dot(bprime);
	    		
	    		c.E(axis0Prime);
	    		c.XE(bprime);
	    		double cComponent = axisNormPrime.dot(c);
	    		
	    		axisNorm.Ea1Tv1(bComponent, b);
	    		axisNorm.PEa1Tv1(cComponent, c);
	    	}
	    		
	    	double dot = axisNorm.dot(siteOrientation[1]);
	    	
	    	if (dot > 0.999999){
	    		u[j+2] = 0;
	    		//return u;
	    	} else
	    	
	    	if (dot < -0.99999){
	    		u[j+2] = Math.PI;
	    		//return u;
	    	} else {
	    	
	    	u[j+2] = Math.acos(dot);
	    	axisNorm.XE(siteOrientation[1]);
	    	
	    	if (axisNorm.dot(siteOrientation[0]) < 0){
	    		u[j+2] = - u[5];
	    		}
	    	}
	    	
	    	j += coordinateDim;
        }
        return u;
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(AtomSet molecules) {
       
    	for (int i=0; i < molecules.getAtomCount() ; i++){
    		
    		IVector3D[] orientation = new IVector3D[3];
    		
    		orientation[0] = (IVector3D)phase.getSpace().makeVector();
    		orientation[1] = (IVector3D)phase.getSpace().makeVector();
    		orientation[2] = (IVector3D)phase.getSpace().makeVector();
    		IAtomGroup molecule = (IAtomGroup)molecules.getAtom(i);
    		
    	    	/*
    	    	 * Determine the Orientation of Each Molecule
    	    	 */
    	    	
    	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().get(0)).getPosition();
    	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().get(5)).getPosition();
    	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().get(10)).getPosition();
    	    	
    	    	
    	    	/*
    	    	 * Determine axis 1 by using Vector Projection
    	    	 */
    	    	axes[0].Ev1Mv2(leafPos5, leafPos0);
    	    	
    	    	axes[1].Ev1Mv2(leafPos10, leafPos0);
    	    	proj.E(axes[0]);
    	    	double dotProd = proj.dot(axes[1]);
    	    	proj.Ea1Tv1(dotProd/axes[0].squared(), axes[0]);
    	    	axes[1].Ev1Mv2(axes[1], proj);
    	    	
    	    	/*
    	    	 * Normalize all the three axes
    	    	 */
    	    	axes[0].TE(1.0/Math.sqrt(axes[0].squared()));
    	    	axes[1].TE(1.0/Math.sqrt(axes[1].squared()));
    	    	axes[2].E(axes[0]);
    	    	axes[2].XE(axes[1]);
    	    	
    	    	orientation[0].E(axes[0]);
    	    	orientation[1].E(axes[1]);
    	    	orientation[2].E(axes[2]);
    	    	
    	    	orientationManager.setAgent(molecule, orientation);
    	    	
    		}
 
    }

    public void setToU(AtomSet molecules, double[] newU) {
        super.setToU(molecules, newU);
        
        int j=3;
        
        for (int i=0; i < molecules.getAtomCount() ; i++){
        	
        	 IAtomGroup molecule = (IAtomGroup)molecules.getAtom(i);
             IVector3D [] siteOrientation = (IVector3D [])orientationManager.getAgent(molecule);
	    	
	    	/*
	    	 *  Finding the component for axes[0]
	    	 *  x = sqrt(1 - u0^2 - u1^2)
	    	 */
	    	
	    	double x = Math.sqrt(1 - u[j]*u[j] - u[j+1]*u[j+1] );

	    	axis0Prime.Ea1Tv1(x, siteOrientation[0]);
	    	axis0Prime.PEa1Tv1(u[j], siteOrientation[1]);
	    	axis0Prime.PEa1Tv1(u[j+1], siteOrientation[2]);
	    	
        	
	    	/*
	    	 * Finding the normal vector of axes[0]
	    	 */
	    	
	    	axisNorm.Ea1Tv1(Math.cos(u[j+2]),siteOrientation[1]);
	    	axisNorm.PEa1Tv1(Math.sin(u[j+2]), siteOrientation[2]);
	    	
	    	/*
	    	 * Treating the rotation angle along axes[0]
	    	 */
	    	
	    	double cosAlpha = axis0Prime.dot(siteOrientation[0]);
	    	
	    	if (cosAlpha > 0.99999999){
	    		axisNormPrime.E(axisNorm);
	    		
	    	} else{
	    		
	    		bprime.Ea1Tv1( 1.0 /cosAlpha, siteOrientation[0]);
	    		bprime.ME(axis0Prime);
	    		bprime.TE(1.0 /Math.sqrt(bprime.squared()));
	    		
	    		b.Ea1Tv1(cosAlpha, bprime);
	    		b.PEa1Tv1(-Math.sqrt(1-cosAlpha*cosAlpha), axis0Prime);
	    		
	    		double bComponent = axisNorm.dot(b);
	    		
	    		c.E(siteOrientation[0]);
	    		c.XE(b);
	    		double cComponent = axisNorm.dot(c);
	    		
	    		axisNormPrime.Ea1Tv1(bComponent, bprime);
	    		axisNormPrime.PEa1Tv1(cComponent, c);
	    	}
	      	
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().get(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().get(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().get(10)).getPosition();
	    	
	    	vector1.Ev1Mv2(leafPos5, leafPos0);
	    	
	     	vector2.Ev1Mv2(leafPos10, leafPos0);
	    	proj.E(vector1);
	    	double dotProd = proj.dot(vector2);
	    	proj.Ea1Tv1(dotProd/vector1.squared(), vector1);
	    	vector2.ME(proj);
	    	
	    	vector1.TE(1/Math.sqrt(vector1.squared()));
	    	vector2.TE(1/Math.sqrt(vector2.squared()));
	    	
	    	/*
	    	 * Comparing the angle between vector1 and axis0Prime 
	    	 */
	    	Normal.E(vector1);
	    	Normal.XE(vector2);
	    	double angle = Math.acos(vector1.dot(axis0Prime));
	    	
	    	
	    	
	    	j += coordinateDim;
        }
       
    }

    private static final long serialVersionUID = 1L;
    protected final IVector3D [] axes;
    protected double [][] basisOrientation ;
    protected final IVector3D com, temp, axis0, axis0Prime; 
    protected final IVector3D proj, proja, projb;
    protected final IVector3D axisNorm, axisNormPrime, b, bprime, c;
    protected final IVector3D vector1, vector2, Normal;
    protected AtomAgentManager orientationManager; 
    protected final AtomGroupAction atomGroupAction;
}