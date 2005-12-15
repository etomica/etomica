/*
 * Created on Jul 12, 2005
 */
package etomica.eam;

import etomica.units.ElectronVolt;

/**
 * Provides the MEAM/EAM parameters for Sn as given in Ravelo & Baskes (1997)
 * and for Ag and Cu as given in Baskes (1992).
 * 
 * Created by K.R. Schadel and D.A. Kofke July 2005.
 */
public class ParameterSetEAM {
	
	public ParameterSetEAM(double E, double A, double beta, double Ro, double alpha,
			double Z, double Zd) {
		this.E = E;
		this.A = A;
		this.beta = beta;
		this.Ro = Ro;
		this.alpha = alpha;
		this.Z = Z;
		this.Zd = Zd;
	}
	
	public final double E;  //cohesive energy (sublimation energy) per atom of reference crystal structure (eV)
	public final double A; //scaling factor for the embedding energy (unitless)
	public final double beta; //exponential decay factor for atomic densities (unitless)
	public final double Ro; //equilibrium nearest-neighbor distance (Angstroms)
	public final double alpha; //exponential decay factor for the universal energy function (unitless)
	public final double Z; //coordination number for the reference crystal stucture (unitless)
	public final double Zd; //coordination number for the true crystal structure (unitless)
	
	public static final ParameterSetEAM Sn = new ParameterSetEAM(ElectronVolt.UNIT.toSim(3.08), 1.0, 6.2, 3.44, 6.20, 12, 6);
	public static final ParameterSetEAM Ag = new ParameterSetEAM(ElectronVolt.UNIT.toSim(2.850), 1.06, 4.46, 2.88, 5.89, 12, 12);
	public static final ParameterSetEAM Cu = new ParameterSetEAM(ElectronVolt.UNIT.toSim(3.540), 1.07, 3.63, 2.56, 5.11, 12, 12);

}
