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
    
    //under revision--- does not work for multiatomics, since intramolecular energy is not considered
    public void thisTrial(Phase phase) {
        double uOld, uNew;
        if(phase.atomCount==0) {return;}
        int i = (int)(rand.nextDouble()*phase.atomCount);
        Atom a = phase.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  //get ith atom in list
        uOld = phase.potentialEnergy.currentValue(a);
        a.displaceWithin(stepSize);
        phase.boundary().centralImage(a.coordinate.position());  //maybe a better way than this
        phase.iterator.moveNotify(a);
        uNew = phase.potentialEnergy.currentValue(a);
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < rand.nextDouble()) {
             a.replace();
             phase.iterator.moveNotify(a);
             return;
        }
        nAccept++;   //accept
    }
}