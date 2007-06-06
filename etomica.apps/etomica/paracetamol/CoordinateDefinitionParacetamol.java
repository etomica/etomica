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
import etomica.atom.AtomAgentManager.AgentSource;
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
        
        x = (IVector3D)phase.getSpace().makeVector();
        y = (IVector3D)phase.getSpace().makeVector();
        z = (IVector3D)phase.getSpace().makeVector();
        xPrime = (IVector3D)phase.getSpace().makeVector();
        zPrime = (IVector3D)phase.getSpace().makeVector();
        arbiVector = (IVector3D)phase.getSpace().makeVector();
        arbiVectorPrime = (IVector3D)phase.getSpace().makeVector();
        
        xNorm = (IVector3D)phase.getSpace().makeVector();
        yNorm = (IVector3D)phase.getSpace().makeVector();
        zNorm = (IVector3D)phase.getSpace().makeVector();
        
        yDoublePrime = (IVector3D)phase.getSpace().makeVector();
        zDoublePrime = (IVector3D)phase.getSpace().makeVector();
        
        xTriplePrime = (IVector3D)phase.getSpace().makeVector();
        yTriplePrime = (IVector3D)phase.getSpace().makeVector();
        zTriplePrime = (IVector3D)phase.getSpace().makeVector();
        
        
        orientationManager = new AtomAgentManager(new OrientationAgentSource(), phase);
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
            }
        	
        	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(t);
            atomGroupAction.actionPerformed(a);
            
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
    
    public void setBasisOrthorhombic(){
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
    	
    }
    
    public void setBasisMonoclinic(){
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
	    	
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
	    	
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
	    	 * 
	    	 */
	      	
	    	u[j] = axis0Prime.dot(siteOrientation[1]);   //components of siteOrientation[1]
	    	u[j+1] = axis0Prime.dot(siteOrientation[2]); //components of siteOrientation[2]
	      	
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
	    	} else
	    	
	    	if (dot < -0.99999){
	    		u[j+2] = Math.PI;
	    	} else {
	    	
	    	u[j+2] = Math.acos(dot);
	    	axisNorm.XE(siteOrientation[1]);
	    	
	    	if (axisNorm.dot(siteOrientation[0]) < 0){
	    		u[j+2] = - Math.acos(dot);  // u[j+2] = -u[5]?????
	    		}
	    	}
	    	
	    	j += coordinateDim/molecules.getAtomCount();
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
    	    	
    	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
    	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
    	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
    	    	
    	    	
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
	    	 *   STEP 1
	    	 * 
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
	    	
	    	/*
	    	 * Determine axis 1 by using Vector Projection
	    	 */
	    	xPrime.Ev1Mv2(leafPos5, leafPos0);
	    	
	    	axisNormPrime.Ev1Mv2(leafPos10, leafPos0);
	    	proj.E(xPrime);
	    	double dotProd = proj.dot(axisNormPrime);
	    	proj.Ea1Tv1(dotProd/xPrime.squared(), xPrime);
	    	axisNormPrime.ME(proj);
	    	
	      	xPrime.TE(1 /Math.sqrt(xPrime.squared()));
	      	axisNormPrime.TE(1 /Math.sqrt(axisNormPrime.squared()));
            
	      	/*
             * Getting the axes for vectorPrime (xPrime, y, zPrime) 
             * 
             * y is the arbitrary rotation axis
             */
	      	
            x.E(siteOrientation[0]);
            y.E(xPrime);
            y.XE(x);
	      	
	      	zPrime.E(xPrime);
	      	zPrime.XE(y);
	      	
	      	/*
	      	 * finding the tensor that brings the arbiVectorPrime to arbiVector
	      	 */ 
	      	double x1Prime = xPrime.x(0);
	      	double x2Prime = xPrime.x(1);
	      	double x3Prime = xPrime.x(2);
	      	
	      	double z1Prime = zPrime.x(0);
	      	double z2Prime = zPrime.x(1);
	      	double z3Prime = zPrime.x(2);
	      	
	      	double x1 = x.x(0);
	      	double x2 = x.x(1);
	      	double x3 = x.x(2);
	      	double y1 = y.x(0);
	      	double y2 = y.x(1);
	      	double y3 = y.x(2);
	      	double z1 = z.x(0);
	      	double z2 = z.x(1);
	      	double z3 = z.x(2);
	      	
	      	
	      	double L11 = x1Prime*x1 + y1*y1 + z1Prime*z1;
	      	double L12 = x2Prime*x1 + y2*y1 + z2Prime*z1;
	      	double L13 = x3Prime*x1 + y3*y1 + z3Prime*z1;
	      	
	      	double L21 = x1Prime*x2 + y1*y2 + z1Prime*z2;
	      	double L22 = x2Prime*x2 + y2*y2 + z2Prime*z2;
	      	double L23 = x3Prime*x2 + y3*y2 + z3Prime*z2;
	      	
	      	double L31 = x1Prime*x3 + y1*y3 + z1Prime*z3;
	      	double L32 = x2Prime*x3 + y2*y3 + z2Prime*z3;
	      	double L33 = x3Prime*x3 + y3*y3 + z3Prime*z3;
	      	
	    	Tensor rotationL = lattice.getSpace().makeTensor();
	    	
	    	rotationL.setComponent(0, 0, L11);
	    	rotationL.setComponent(0, 1, L12);
	    	rotationL.setComponent(0, 2, L13);
	    	rotationL.setComponent(1, 0, L21);
	    	rotationL.setComponent(1, 1, L22);
	    	rotationL.setComponent(1, 2, L23);
	    	rotationL.setComponent(2, 0, L31);
	    	rotationL.setComponent(2, 1, L32);
	    	rotationL.setComponent(2, 2, L33);
	    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationL);
	        atomGroupAction.actionPerformed(molecule);
	      	
	      	/*
	      	 *  STEP 2
	      	 * 
	      	 * Treat x as the new vector along the x-axis (vector formed from atom0 to atom10)
	      	 * and yNorm is normal to x-axis
	      	 */
	        
	    	xNorm.Ev1Mv2(leafPos5, leafPos0);
	    	
	    	yNorm.Ev1Mv2(leafPos10, leafPos0);
	    	proj.E(xNorm);
	    	double dotProdNorm = proj.dot(yNorm);
	    	proj.Ea1Tv1(dotProdNorm/xNorm.squared(), xNorm);
	    	yNorm.ME(proj);
	    	
	      	xNorm.TE(1 /Math.sqrt(xNorm.squared()));
	      	yNorm.TE(1 /Math.sqrt(yNorm.squared()));
	        
	        zNorm.E(xNorm);
	        zNorm.XE(yNorm);
	        
	        //Getting the DoublePrime axes
	    	
	    	yDoublePrime.Ea1Tv1(Math.cos(u[j+2]),siteOrientation[1]);
	    	yDoublePrime.PEa1Tv1(Math.sin(u[j+2]), siteOrientation[2]);
	        
	        zDoublePrime.E(xNorm);
	        zDoublePrime.XE(yDoublePrime);
	        
	      	/*
	      	 * finding the tensor that brings the arbiVector to arbiVectorDoublePrime
	      	 */ 
	      	double x1Norm = xNorm.x(0);
	      	double x2Norm = xNorm.x(1);
	      	double x3Norm = xNorm.x(2);
	      	
	      	double y1Norm = yNorm.x(0);
	      	double y2Norm = yNorm.x(1);
	      	double y3Norm = yNorm.x(2);
	      	
	      	double z1Norm = zNorm.x(0);
	      	double z2Norm = zNorm.x(1);
	      	double z3Norm = zNorm.x(2);

	      	double y1DoublePrime = yDoublePrime.x(0);
	      	double y2DoublePrime = yDoublePrime.x(1);
	      	double y3DoublePrime = yDoublePrime.x(2);
	      	double z1DoublePrime = zDoublePrime.x(0);
	      	double z2DoublePrime = zDoublePrime.x(1);
	      	double z3DoublePrime = zDoublePrime.x(2);
	      	
	      	
	      	double M11 = x1Norm*x1Norm + y1Norm*y1DoublePrime + z1Norm*z1DoublePrime;
	      	double M12 = x2Norm*x1Norm + y2Norm*y1DoublePrime + z2Norm*z1DoublePrime;
	      	double M13 = x3Norm*x1Norm + y3Norm*y1DoublePrime + z3Norm*z1DoublePrime;
	      	
	      	double M21 = x1Norm*x2Norm + y1Norm*y2DoublePrime + z1Norm*z2DoublePrime;
	      	double M22 = x2Norm*x2Norm + y2Norm*y2DoublePrime + z2Norm*z2DoublePrime;
	      	double M23 = x3Norm*x2Norm + y3Norm*y2DoublePrime + z3Norm*z2DoublePrime;
	      	
	      	double M31 = x1Norm*x3Norm + y1Norm*y3DoublePrime + z1Norm*z3DoublePrime;
	      	double M32 = x2Norm*x3Norm + y2Norm*y3DoublePrime + z2Norm*z3DoublePrime;
	      	double M33 = x3Norm*x3Norm + y3Norm*y3DoublePrime + z3Norm*z3DoublePrime;
	      	
	    	Tensor rotationM = lattice.getSpace().makeTensor();
	    	
	    	rotationL.setComponent(0, 0, M11);
	    	rotationL.setComponent(0, 1, M12);
	    	rotationL.setComponent(0, 2, M13);
	    	rotationL.setComponent(1, 0, M21);
	    	rotationL.setComponent(1, 1, M22);
	    	rotationL.setComponent(1, 2, M23);
	    	rotationL.setComponent(2, 0, M31);
	    	rotationL.setComponent(2, 1, M32);
	    	rotationL.setComponent(2, 2, M33);
	    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationM);
	        atomGroupAction.actionPerformed(molecule);
	        
	    	
	    	/*
	    	 *    STEP  3
	    	 * 
	    	 * 
	    	 * 
	    	 *  First we find the component for xTriplePrime
	    	 *  x = sqrt(1 - u[j]^2 - u[j+1]^2)
	    	 */
	    	
	    	double x = Math.sqrt(1 - u[j]*u[j] - u[j+1]*u[j+1]);

	    	xTriplePrime.Ea1Tv1(x, siteOrientation[0]);
	    	xTriplePrime.PEa1Tv1(u[j], siteOrientation[1]);
	    	xTriplePrime.PEa1Tv1(u[j+1], siteOrientation[2]);
	    	
	    	yTriplePrime.E(xTriplePrime);
	    	yTriplePrime.XE(xNorm);
	    
	    	zTriplePrime.E(xTriplePrime);
	    	zTriplePrime.XE(yTriplePrime);
	    	
	      	/*
	      	 * finding the tensor that brings the arbiVectorDoublePrime to arbiVectorTriplePrime
	      	 */ 
	      	double x1TriplePrime = xTriplePrime.x(0);
	      	double x2TriplePrime = xTriplePrime.x(1);
	      	double x3TriplePrime = xTriplePrime.x(2);
	    	
	      	double y1TriplePrime = yTriplePrime.x(0);
	      	double y2TriplePrime = yTriplePrime.x(1);
	      	double y3TriplePrime = yTriplePrime.x(2);
	      	
	      	double z1TriplePrime = zTriplePrime.x(0);
	      	double z2TriplePrime = zTriplePrime.x(1);
	      	double z3TriplePrime = zTriplePrime.x(2);
	      	
	      	
	      	double N11 = x1Norm*x1TriplePrime + y1DoublePrime*y1TriplePrime + z1DoublePrime*z1TriplePrime;
	      	double N12 = x2Norm*x1TriplePrime + y2DoublePrime*y1TriplePrime + z2DoublePrime*z1TriplePrime;
	      	double N13 = x3Norm*x1TriplePrime + y3DoublePrime*y1TriplePrime + z3DoublePrime*z1TriplePrime;
	      	
	      	double N21 = x1Norm*x2TriplePrime + y1DoublePrime*y2TriplePrime + z1DoublePrime*z2TriplePrime;
	      	double N22 = x2Norm*x2TriplePrime + y2DoublePrime*y2TriplePrime + z2DoublePrime*z2TriplePrime;
	      	double N23 = x3Norm*x2TriplePrime + y3DoublePrime*y2TriplePrime + z3DoublePrime*z2TriplePrime;
	      	
	      	double N31 = x1Norm*x3TriplePrime + y1DoublePrime*y3TriplePrime + z1DoublePrime*z3TriplePrime;
	      	double N32 = x2Norm*x3TriplePrime + y2DoublePrime*y3TriplePrime + z2DoublePrime*z3TriplePrime;
	      	double N33 = x3Norm*x3TriplePrime + y3DoublePrime*y3TriplePrime + z3DoublePrime*z3TriplePrime;
	      	
	    	Tensor rotationN = lattice.getSpace().makeTensor();
	    	
	    	rotationL.setComponent(0, 0, N11);
	    	rotationL.setComponent(0, 1, N12);
	    	rotationL.setComponent(0, 2, N13);
	    	rotationL.setComponent(1, 0, N21);
	    	rotationL.setComponent(1, 1, N22);
	    	rotationL.setComponent(1, 2, N23);
	    	rotationL.setComponent(2, 0, N31);
	    	rotationL.setComponent(2, 1, N32);
	    	rotationL.setComponent(2, 2, N33);
	    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationN);
	        atomGroupAction.actionPerformed(molecule);

	    	
	    	j += coordinateDim/molecules.getAtomCount();
        }
       
    }

    private static final long serialVersionUID = 1L;
    protected final IVector3D [] axes;
    protected double [][] basisOrientation ;
    protected final IVector3D com, temp, axis0, axis0Prime; 
    protected final IVector3D proj, proja, projb;
    protected final IVector3D axisNorm, axisNormPrime, b, bprime, c;
    protected final IVector3D x, y, z, xPrime, zPrime;
    protected final IVector3D arbiVector, arbiVectorPrime;
    protected final IVector3D xNorm, yNorm, zNorm;
    protected final IVector3D yDoublePrime, zDoublePrime;
    protected final IVector3D xTriplePrime, yTriplePrime, zTriplePrime;
    
    
    protected AtomAgentManager orientationManager; 
    protected final AtomGroupAction atomGroupAction;

    protected static class OrientationAgentSource implements AgentSource, Serializable {
        
        public OrientationAgentSource() {
        }
        public Class getAgentClass() {
            return IVector [].class;
        }
        public Object makeAgent(IAtom atom) {
            return null;
        }
        public void releaseAgent(Object agent, IAtom atom) {
            //nothing to do
        }
        
        private static final long serialVersionUID = 1L;
    }
    
}
