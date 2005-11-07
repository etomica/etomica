package etomica.modules.dcvgcmd;

import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.units.Dimension;

public class MeterFlux extends DataSourceScalar implements Meter {
    
	public MeterFlux(MyMCMove move, IntegratorDCVGCMD integrator) {
		super("Flux",Dimension.UNDEFINED);
		mcMove = move;
		integratorMD = integrator;
	}
 
	public double getDataAsScalar() {
		double t1 = integratorMD.elapsedTime();
        if(t1 == t0) return Double.NaN;
		int n1 = mcMove.getDeltaN();
		double rate = (n1 - n0)/(t1 - t0)/(phase.getBoundary().dimensions().x(0)*phase.getBoundary().dimensions().x(1));
		n0 = n1;
		t0 = t1;
		return rate;
	}
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private double t0 = 0;
    private int n0 = 0;
    private MyMCMove mcMove;
    private IntegratorDCVGCMD integratorMD;
}