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
	public MeterVirialB3(Simulation sim, double refSigma, P0Cluster simulationPotential, Potential2 targetPotential) {
		super(sim,
				1.0, //refTemperature
				new etomica.virial.cluster.C3(new MayerHardSphere(refSigma)),//refCluster 
				B3HS(refSigma),//refIntegral
				B3Clusters(targetPotential),//target clusters
				simulationPotential);//sampling potential
		System.out.println("In MeterVirialB3, refIntegral = "+B3HS(refSigma));
	}
	
	public static double B3HS(double sigma) {
		return etomica.virial.cluster.Standard.C3HS(sigma);
	}
	
	public static Cluster[] B3Clusters(Potential2 potential) {
		return new Cluster[] {new Ring(3, -1./3., new MayerGeneral(potential))};
	}

}
