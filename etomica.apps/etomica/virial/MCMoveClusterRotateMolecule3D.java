package etomica.virial;

import etomica.action.MoleculeAction;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.space.ISpace;


public class MCMoveClusterRotateMolecule3D extends MCMoveRotateMolecule3D {

    public MCMoveClusterRotateMolecule3D(IPotentialMaster potentialMaster,
            IRandom random, ISpace _space) {
        super(potentialMaster, random, _space);
        weightMeter = new MeterClusterWeight(potential);
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        weightMeter.setBox(p);
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        while (molecule.getIndex() == 0) {
            molecule = moleculeSource.getMolecule();
        }
        uOld = weightMeter.getDataAsScalar();
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(random.nextInt(3),dTheta);

        IAtomList leafAtoms = molecule.getChildList();
        IAtom first = leafAtoms.getAtom(0);
        r0.E(first.getPosition());
        doTransform();

        if (trialCount-- == 0) {
            relaxAction.actionPerformed(molecule);
            trialCount = relaxInterval;
        }

        uNew = weightMeter.getDataAsScalar();
        ((BoxCluster)box).trialNotify();
        return true;
    }
    
    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        return uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((BoxCluster)box).acceptNotify();
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((BoxCluster)box).rejectNotify();
    }
    
    public void setRelaxAction(MoleculeAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    protected int trialCount, relaxInterval = 100;
    protected MoleculeAction relaxAction;
}
