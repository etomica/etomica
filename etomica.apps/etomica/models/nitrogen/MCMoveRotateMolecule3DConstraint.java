package etomica.models.nitrogen;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.ISpace;
import etomica.space.RotationTensor;


/**
 * MCMoveRotate that moves the molecules according to pre-set constraint angle
 *  which is given by the parameter "angle"
 *  
 * This class does not consider rotational energy, so the variable "energyChange"
 *  always returns zero
 *  
 * getB() always returns positive one is to ensure the proposed move is always 
 *  accepted
 * 
 * @author Tai Boon Tan
 *
 */
public class MCMoveRotateMolecule3DConstraint extends MCMoveMolecule {
    
    private static final long serialVersionUID = 2L;
    protected transient IVectorMutable r0;
    protected transient RotationTensor rotationTensor;
    protected IAtomPositionDefinition positionDefinition;
    protected double constraintAngle;
    
    public MCMoveRotateMolecule3DConstraint(IPotentialMaster potentialMaster, IRandom random,
    		                      ISpace _space, double angle) {
        super(potentialMaster, random, _space, 0.5*Math.PI, Math.PI);
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        positionDefinition = new AtomPositionGeometricCenter(space);
        constraintAngle = angle;
    }
     
    public boolean doTrial() {
//        System.out.println("doTrial MCMoveRotateMolecule called");
        
        if(box.getMoleculeList().getMoleculeCount()==0) {molecule = null; return false;}
            
        molecule = moleculeSource.getMolecule();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        double randomValue = (2*random.nextDouble() - 1.0);
        double constant = (constraintAngle/180.0) *Math.PI;
        double dTheta = randomValue*constant;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

        r0.E(positionDefinition.position(molecule));
        doTransform();
        
        energyMeter.setTarget(molecule);
        return true;
    }
    
    protected void doTransform() {
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
        rotationTensor.invert();
        doTransform();
    }
    
    public double getB() {
//    	double energy = energyMeter.getDataAsScalar();
//    	if(Double.isInfinite(energy)){
//    		return -1.0;
//    	}
        return 1.0; //always accept the rotational move, See the acceptance criteria in IntegratorMC
    }
    
    public double energyChange() { return 0.0;}
}
