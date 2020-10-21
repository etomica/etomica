package etomica.spin.heisenberg;

import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveRotatePair extends MCMoveRotate{
    protected final IPotentialAtomic p2;
    public MCMoveRotatePair(IPotentialAtomic p2, IRandom random, Space _space) {
        super(null, random, _space);
        this.p2=p2;

    }

    public boolean doTrial() {
        if (box.getMoleculeList().size() == 0) {
            return false;
        }
        atom = atomSource.getAtom();
//        atom = box.getLeafList().getAtom(0);//TODO only rotate first molecule!!

        energyMeter.setTarget(atom);
        uOld = p2.energy(box.getLeafList());
        iOrientation = ((IAtomOriented)atom).getOrientation();
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);
        /*if (atom.getLeafIndex()==0){
            iOrientation.getDirection().E(new double[]{-1,0});
        }
        else {
            iOrientation.getDirection().E(new double[]{1,0});

        }*/
        return true;
    }
    public double getChi(double temperature) {
        uNew = p2.energy(box.getLeafList());
        double chi = Math.exp(-(uNew - uOld) / temperature);
//        System.out.println("chi "+chi+" "+uOld+" "+uNew+" "+temperature);
        return chi;
    }
}
