package simulate;

import java.util.Random;

public class MCMoveVolume extends MCMove {
    
    private transient final double[] dr = new double[Space.D];
    private final Random rand = new Random();
    double pressure;

    public MCMoveVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(0.0);
    }
    
    public void thisTrial(Phase phase) {
        double hOld, hNew, vOld, vNew;
        vOld = phase.space.volume;
        hOld = phase.potentialEnergy() + pressure*vOld*Constants.PV2T;
        double vScale = (2.*rand.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)Space.D);
        phase.space.inflate(rScale);
        for(Molecule m=phase.firstMolecule; m!=null; m=m.getNextMolecule()) {
            m.inflate(rScale);
        }
        hNew = phase.potentialEnergy() + pressure*vNew*Constants.PV2T;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.nMoleculeTotal+1)*vScale)
                < rand.nextDouble()) 
            {  //reject
              phase.space.inflate(1.0/rScale);
              for(Molecule m=phase.firstMolecule; m!=null; m=m.getNextMolecule()) {m.replace();}
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public double getPressure() {return pressure;}
    public final void setLogPressure(int lp) {pressure = Math.pow(10.,(double)lp);}
}