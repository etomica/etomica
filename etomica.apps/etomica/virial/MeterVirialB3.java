package etomica.virial;

import etomica.*;
import etomica.virial.cluster.Ring;

/**
 * @author David Kofke
 *
 */
public class MeterVirialB3 extends MeterVirial {

	/**
	 * Constructor for MeterVirialB3.
	 * @param sim
	 * @param pairSet
	 * @param sigma
	 * @param simulationPotential
	 * @param targetPotential
	 */
	public MeterVirialB3(Simulation sim, PairSet pairSet, double refSigma, P2Cluster simulationPotential, Potential2 targetPotential) {
		super(sim, pairSet, 
				1.0, new etomica.virial.cluster.C3(new MayerHardSphere(refSigma)), B3HS(refSigma),
				B3Clusters(targetPotential),
				simulationPotential);
	}
	
	public static double B3HS(double sigma) {
		double b0 = 2.*Math.PI/3 * sigma*sigma*sigma;
		return -5./8. * b0 * b0;
	}
	
	public static Cluster[] B3Clusters(Potential2 potential) {
		return new Cluster[] {new Ring(3, -1./3., new MayerGeneral(potential))};
	}

}
