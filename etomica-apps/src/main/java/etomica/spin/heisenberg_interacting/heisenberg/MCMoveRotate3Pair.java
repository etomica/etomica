package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomPair;
import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveRotate3Pair extends MCMoveRotate {
    protected final IPotentialAtomic p2;

    public MCMoveRotate3Pair(IPotentialAtomic p2, IRandom random, Space _space) {
        super(null, random, _space);
        this.p2 = p2;

    }

    public boolean doTrial() {
        if (box.getMoleculeList().getMoleculeCount() == 0) {
            return false;
        }
        atom = atomSource.getAtom();
        energyMeter.setTarget(atom);

        AtomPair pair = new AtomPair();
        pair.atom0 = box.getLeafList().getAtom(0);
        pair.atom1 = box.getLeafList().getAtom(1);//01
        uOld = p2.energy(pair);
        pair.atom1 = box.getLeafList().getAtom(2);//02
        uOld += p2.energy(pair);
        pair.atom0 = box.getLeafList().getAtom(1);//12
        uOld += p2.energy(pair);

        iOrientation = ((IAtomOriented) atom).getOrientation();
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
        AtomPair pair = new AtomPair();
        pair.atom0 = box.getLeafList().getAtom(0);
        pair.atom1 = box.getLeafList().getAtom(1);//01
        uNew = p2.energy(pair);
        pair.atom1 = box.getLeafList().getAtom(2);//02
        uNew += p2.energy(pair);
        pair.atom0 = box.getLeafList().getAtom(1);//12
        uNew += p2.energy(pair);
        double chi = Math.exp(-(uNew - uOld) / temperature);
//        System.out.println("chi "+chi+" "+uOld+" "+uNew+" "+temperature);
        return chi;
    }
}
