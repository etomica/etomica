package etomica.modules.dcvgcmd;

import etomica.EtomicaElement;
import etomica.Phase;
import etomica.data.meter.MeterScalar;
import etomica.units.Dimension;

public class MeterFlux extends MeterScalar implements EtomicaElement{
    
	double t0 = 0;
	int n0 = 0;
	MyMCMove mcMove;
	IntegratorDCVGCMD integratorMD;
	
	
	//private final MeterKineticEnergy meterKE;
    
	public MeterFlux(MyMCMove move, IntegratorDCVGCMD integrator) {
		super();
		setLabel("Flux");
		mcMove = move;
		integratorMD = integrator;
	}
 
	public double getDataAsScalar(Phase phase) {
		double t1 = integratorMD.elapsedTime();
		int n1 = mcMove.getDeltaN();
		double rate = (n1 - n0)/(t1 - t0)/(phase.boundary().dimensions().x(0)*phase.boundary().dimensions().x(1));
		n0 = n1;
		t0 = t1;
		return rate;
	}
    
	public Dimension getDimension() {
		return Dimension.UNDEFINED;}

}