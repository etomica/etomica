package etomica.chem.models.water;
import etomica.units.*;

/**
 * @author kofke
 * The simple-point-charge model of water.
 */
public class SPC extends Abstract3Site {

	/**
	 * Constructor for SPC.
	 * sigma = 3.1670 Angstroms; epsilon/kB = 78.23 K; qH = 0.41 e; r0H = 1.0
	 * Angstroms; thetaHOH = 109.5 degrees
	 */
	public SPC() {
		super(3.1670, Kelvin.UNIT.toSim(78.23), 0.41, 1.0, Degree.UNIT.toSim(109.5));//sigma, epsilon, qH, rOH, theta
	}

}
