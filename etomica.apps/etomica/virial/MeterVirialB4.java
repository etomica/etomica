package etomica.virial;

import etomica.*;
import etomica.potential.Potential2;
import etomica.virial.cluster.ReeHoover;
import etomica.virial.cluster.Standard;

/**
 * @author kofke
 *
 * Meter to compute the 4th virial coefficient.
 * allMeters[0] is the sum of the cluster diagrams, while allMeters[1], [2], [3]
 * correspond to the three cluster diagrams.
 * 
 */
public class MeterVirialB4 extends MeterVirial {
	
	public MeterVirialB4(Simulation sim, double refSigma, P0Cluster simulationPotential, Potential2 targetPotential) {
		super(sim, 
				1.0, new etomica.virial.cluster.D4(new MayerHardSphere(refSigma)), D4HS(refSigma),
				B4Clusters(targetPotential),
				simulationPotential);
	}
	
	public static double D4HS(double sigma) {
		double b0 = 2.*Math.PI/3 * sigma*sigma*sigma;
		return -0.9714 * b0 * b0 * b0;
	}
	
	public static Cluster[] B4Clusters(Potential2 potential) {
		MayerFunction f = new MayerGeneral(potential);
		Cluster c1 = new ReeHoover(4, 3.0, new Cluster.BondGroup(f, Standard.D4));
		Cluster c2 = new Cluster(4, -2.0, new Cluster.BondGroup(f, Standard.D6));
		return new Cluster[] {c1, c2};
	}
	
}
