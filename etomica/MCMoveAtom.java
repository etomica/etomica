package simulate;

import java.util.Random;

public class MCMoveAtom extends MCMove {
    
    private final Random rand = new Random();

    public MCMoveAtom() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
    }
    
    public void thisTrial(PhaseSpace phaseSpace) {
        double uOld, uNew;
        if(phaseSpace.atomCount==0) {return;}
        int i = (int)(rand.nextDouble()*phaseSpace.nAtomTotal);
        Atom a = phaseSpace.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  //get ith atom in list
        uOld = phaseSpace.potentialEnergy.currentValue(a);
        a.coordinate.displaceWithin(stepSize);
        uNew = phaseSpace.potentialEnergy.currentValue(a);
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < rand.nextDouble()) {
             a.coordinate.replace();
             return;
        }
        nAccept++;   //accept
    }
}