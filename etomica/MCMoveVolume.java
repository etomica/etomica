package etomica;

import java.util.Random;
import etomica.units.Dimension;

public class MCMoveVolume extends MCMove {
    
    private final Random rand = new Random();
    protected double pressure;
    protected PhaseAction.Inflate inflate;

    public MCMoveVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        inflate = new PhaseAction.Inflate(phase);
    }
    
    public void thisTrial() {
        double hOld, hNew, vOld, vNew;
        vOld = phase.volume();
        hOld = phase.energy.potential() + pressure*vOld;
        double vScale = (2.*Math.random()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)phase.parentSimulation().space().D());
        inflate.actionPerformed(phase,rScale);
        hNew = phase.energy.potential() + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.moleculeCount+1)*vScale)
                < Math.random()) 
            {  //reject
              inflate.retractAction();
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
}