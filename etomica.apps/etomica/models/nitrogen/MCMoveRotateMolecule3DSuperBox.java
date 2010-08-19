package etomica.models.nitrogen;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.ISpace;
import etomica.space.RotationTensor;


/**
 * MC Rotate Move 3D for superbox
 * 
 * 
 * @author taitan
 *
 */
public class MCMoveRotateMolecule3DSuperBox extends MCMoveMolecule {
    
    private static final long serialVersionUID = 2L;
    protected transient IVectorMutable r0;
    protected transient RotationTensor rotationTensor;
    protected IAtomPositionDefinition positionDefinition;
    public int count;
    public int count1;
    public boolean flag = false;
    public boolean flag1 = false;
    protected int nA;
    protected int[][] molIndex;
    protected int molNum;
    protected IRandom random;
    protected BasisCell[] basisCell;
    protected IMolecule affectedMol;
    protected CoordinateDefinitionNitrogenSuperBox coordinateDef;
    
    public MCMoveRotateMolecule3DSuperBox(IPotentialMaster potentialMaster, IRandom random,
    		                      ISpace _space, int nC, int basis, CoordinateDefinitionNitrogenSuperBox coordinateDef) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI);
        this.basisCell = coordinateDef.getBasisCells();
        this.random = random;
        this.coordinateDef = coordinateDef;
        
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        positionDefinition = new AtomPositionGeometricCenter(space);
        
        nA = (nC*nC*nC)*basis;
        molIndex = new int[27][nA];
        
        int layerNum = (nC*basis*3)*(3*nC);
        int oneThirdLayerNum = (nC*basis*3)*(nC);
        
        int axisNum = nC*basis*3;
        int cellAxisNum = nC*basis;
        int ix= 0;
        int iy= 0;
        int iz= 0;
        
        for (int iCell=0; iCell<molIndex.length; iCell++){
        	
        	int counter = 0;
        	int cellConst = 0;
        	
        	if(iCell%27 > 17){
        		ix = 2;
        	} else if (iCell%27 > 8 && iCell%27 <=17){
        		ix = 1;
        	} else {
        		ix = 0;
        	}
        	
        	if(iCell%9 > 5){
        		iy = 2;
        	} else if (iCell%9 > 2 && iCell%9 <=5){
        		iy = 1;
        	} else {
        		iy = 0;
        	}
        	
        	iz = iCell%3;
        
        	cellConst = iz*(cellAxisNum) + iy*oneThirdLayerNum + ix*nC*layerNum;
        	
	        for(int xnC=0; xnC<nC; xnC++){
		        for(int ynC=0; ynC<nC; ynC++){
			        for(int i=0; i<nC*basis; i++){
			        	molIndex[iCell][counter] = (i+(ynC*axisNum)+(xnC*layerNum)+cellConst);
			        	++counter;
			        }
		        }
	        }
        }
        
        
    }
     
    public boolean doTrial() {
//        System.out.println("doTrial MCMoveRotateMolecule called");
        
        if(box.getMoleculeList().getMoleculeCount()==0) {molecule = null; return false;}
            
        molNum = random.nextInt(nA);
        molecule = basisCell[0].molecules.getMolecule(molNum);
        
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

        for(int i=0; i<molIndex.length; i++){
        	affectedMol = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
	        r0.E(coordinateDef.getLatticePosition(affectedMol));
	        doTransform(affectedMol, r0);
        }
        
        energyMeter.setTarget(molecule);
        return true;
    }
    
    protected void doTransform(IMolecule molecule, IVector r0) {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            IVectorMutable r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }
    
    public void rejectNotify() {
    	
        for(int i=0; i<molIndex.length; i++){
        	affectedMol = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
        	r0.E(coordinateDef.getLatticePosition(affectedMol));
	        rotationTensor.invert();
	        doTransform(affectedMol, r0);
        }
        
    }
}
