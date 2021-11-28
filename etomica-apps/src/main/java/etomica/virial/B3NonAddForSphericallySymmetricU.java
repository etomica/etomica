/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.potential.IPotential2;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.potential.Potential3Soft;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;


/**
 * 
 * Computes the classical, non-additive component of B3 for a spherically symmetric potential using quadrature.
 * 
 * @author kate
 *
 */



public class B3NonAddForSphericallySymmetricU {
	
	public B3NonAddForSphericallySymmetricU() {
		
	}
	
public static void main(String[] args) {

	
	Space space = Space3D.getInstance();
	
	double temperature = 25; // Kelvin
	if (args.length == 6 ) {

		 r02Min = Double.parseDouble(args[0]); //Bohr radii
		 r02Max = Double.parseDouble(args[1]); //Bohr radii
		 delr02 = Double.parseDouble(args[2]); //Bohr radii
		 r03Max = Double.parseDouble(args[3]); //Bohr radii
		 delcostheta=Double.parseDouble(args[4]); //Radians
		 temperature=Double.parseDouble(args[5]); //Kelvin
		 
    } 

	
	double[] temps = new double[] { temperature }; // Kelvin
	
	P2HePCKLJS p2 = new P2HePCKLJS();
	
	P3CPSNonAdditiveHe p3 = new P3CPSNonAdditiveHe(space);
	
	System.out.println("T(K)    B3NonAdd (cm6/mol2)");
	
	
	System.out.println();
	
	for (int t=0; t<temps.length; t++) {
		
		double temp = temps[t];
  
		 double B3NonAdd = computeB3NonAdd(p2,p3,space, temp);
			
		System.out.println(temp + "    "  +  B3NonAdd);

	}
	

	}

	public static double computeB3NonAdd(IPotential2 p2, Potential3Soft p3, Space space, double temp) {

		Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        Atom atom3 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
	
        double costhetaMax = 1.0;
        
        
        double r03Min = 0; 
        double costhetaMin = 0;
        
        
        double delr03 = delr02; //Bohr radii
        
        
        long Nr02 = Math.round((r02Max-r02Min)/delr02);
        long Nr03 = Math.round((r03Max-r03Min)/delr03);
        
        
        r02Min = r02Min*AngstromPerBohrRadius;
		r02Max = r02Max*AngstromPerBohrRadius;
		delr02 = delr02*AngstromPerBohrRadius;
		r03Max = r03Max*AngstromPerBohrRadius;
        r03Min = r03Min*AngstromPerBohrRadius;
        
        long Ncostheta = Math.round((costhetaMax-costhetaMin)/delcostheta);

		double B3NonAdd = 0;
		System.out.println("r12Max (a0)    B3NonAdd (cm6/mol2)");

		for (int nr02=0; nr02 <= Nr02; nr02++) {
			
			double r02= r02Min + nr02*delr02; 
			double r12 = 2.0*r02;
			
			double Icostheta = 0;
			for (int ncostheta=0; ncostheta <= Ncostheta; ncostheta++) {
				
				// angle between origin and position of atom 3
				double costheta = costhetaMin + delcostheta*ncostheta;
			
				double Ir03 = 0;
				for (int nr03=0; nr03 <= Nr03; nr03++) {
					
					double r03 = r03Min + nr03*delr03;

					// origin halfway between atom 1 and atom 2
					double x13 = r03*costheta;
					double y13 = r03*Math.sqrt(1.0-costheta*costheta);
					
					double integrand = 0;
								
					Vector r1 = Vector.of(new double[]{-r02, 0, 0});
			        Vector r2 = Vector.of(new double[]{r02, 0, 0});
			        Vector r3 = Vector.of(new double[]{x13, y13, 0});
					double r212 = 4*r02*r02;
					double r213 = (x13+r02)*(x13+r02) + y13*y13;
					double r223 = (x13-r02)*(x13-r02) + y13*y13;
			        
			        atom1.getPosition().E(r1);
			        atom2.getPosition().E(r2);
			        atom3.getPosition().E(r3);
			        
			        
			        double r13 = Math.sqrt((r1.getX(0)-r3.getX(0))*(r1.getX(0)-r3.getX(0))+(r1.getX(1)-r3.getX(1))*(r1.getX(1)-r3.getX(1)));
			        double r23 = Math.sqrt((r2.getX(0)-r3.getX(0))*(r2.getX(0)-r3.getX(0))+(r2.getX(1)-r3.getX(1))*(r2.getX(1)-r3.getX(1)));
					double u123Add = p2.u(r12*r12) + p2.u(r13*r13) + p2.u(r23*r23); //Kelvin
					double u123NonAdd = p3.u(r212, r213, r223);
					double e123 = Math.exp(-(u123NonAdd+u123Add)/temp);
					double e123Add = Math.exp(-u123Add/temp);			

					integrand = e123-e123Add;
					
					if (Double.isNaN(integrand)) {
						integrand = e123-e123Add;
					}
					
					
					if (nr03 == 0 | nr03 == Nr03) {
						
						Ir03 = Ir03 + 0.5*integrand*(2*Math.PI*r03)*r03*delr03;
						
					} else {
						
						Ir03 = Ir03 + integrand*(2*Math.PI*r03)*r03*delr03;
						
					}
								
				}
				
				if (ncostheta == 0 | ncostheta == Ncostheta) {
				
					
					Icostheta = Icostheta + 0.5*Ir03*2.0*delcostheta;
			
				} else {
					
					Icostheta = Icostheta +     Ir03*2.0*delcostheta;
					
				}

			}
			
			if (nr02 == 0 | nr02 == Nr02) {
				B3NonAdd = B3NonAdd + 0.5*Icostheta*(4.0*Math.PI*r12*r12*2.0*delr02)/(-3.0);
			} else {
				B3NonAdd = B3NonAdd +     Icostheta*(4.0*Math.PI*r12*r12*2.0*delr02)/(-3.0);
			}
			
			System.out.println(r02/AngstromPerBohrRadius+" " + B3NonAdd*0.60221415*0.60221415);
		}

		return B3NonAdd*0.60221415*0.60221415; //cm6/mol2
	 }
	
	 static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
	 static double delcostheta=0.005;
	 static double delr02 = 0.05; //Bohr radii
	 static double r02Min = 1.0;  //Bohr radii
	 static double r02Max = 20;   //Bohr radii
	 static double r03Max = 20;   //Bohr radii
}
