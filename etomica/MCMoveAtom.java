package etomica;

import java.util.Random;
import etomica.units.Dimension;

public class MCMoveAtom extends MCMove {
    
    private final Random rand = new Random();

    public MCMoveAtom() {
        super();
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    
    
    //under revision--- does not work for multiatomics, since intramolecular energy is not considered
    public void thisTrial() {
        double uOld, uNew;
        if(phase.atomCount==0) {return;}
        int i = (int)(rand.nextDouble()*phase.atomCount);
        Atom a = phase.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  //get ith atom in list
        uOld = phase.potential.energy(a);
        a.displaceWithin(stepSize);
        phase.boundary().centralImage(a.coordinate.position());  //maybe a better way than this
        phase.iteratorFactory().moveNotify(a);
        uNew = phase.potential.energy(a);
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < rand.nextDouble()) {
             a.replace();
             phase.iteratorFactory().moveNotify(a);
             return;
        }
        nAccept++;   //accept
    }
}