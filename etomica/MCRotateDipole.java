package etomica;

import java.util.Random;
import etomica.electrostatics.*;

//written for 2D system
public class MCRotateDipole extends MCMove {
    
    private final Random rand = new Random();

    public MCRotateDipole() {
        super();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/10.0);
    }
    
    //under revision--- does not work for multiatomics, since intramolecular energy is not considered
    public void thisTrial() {
        double uOld, uNew;
        if(phase.atomCount==0) {return;}
        int i = (int)(rand.nextDouble()*phase.atomCount);
        Atom a = phase.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  //get ith atom in list
        uOld = phase.energy.meterPotential().currentValue(a);
        
        Dipole d = (Dipole)a.type.electroType();
        Space2D.Vector e = (Space2D.Vector)d.e();
        double theta = Math.atan2(e.y,e.x);
        double delta = (2.0*rand.nextDouble()-1.0)*stepSize;
        theta += delta;
        e.x = Math.cos(theta);
        e.y = Math.sin(theta);
        uNew = phase.energy.meterPotential().currentValue(a);
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