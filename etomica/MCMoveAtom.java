package simulate;

import java.util.Random;

public class MCMoveAtom extends MCMove {
    
    private transient final double[] dr = new double[Space.D];
    private final Random rand = new Random();

    public MCMoveAtom() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
    }
    
    public void thisTrial(Phase phase) {
        double uOld, uNew;
        if(phase.nAtomTotal==0) {return;}
        int i = (int)(rand.nextDouble()*phase.nAtomTotal);
        Atom a = phase.firstAtom;
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.getNextAtom();}  //get ith atom in list
        uOld = a.potentialEnergy();
        Space.randomVector(dr, stepSize, rand);
        a.displace(dr);
        uNew = a.potentialEnergy();
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < rand.nextDouble()) {
             a.replace();
             return;
        }
        nAccept++;   //accept
    }
}