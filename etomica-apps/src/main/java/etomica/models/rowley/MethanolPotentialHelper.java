/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.atom.IAtomType;
import etomica.api.IPotentialAtomic;
import etomica.atom.iterator.ApiBuilder;
import etomica.potential.P2ModifiedMorse;
import etomica.potential.P2Morse;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.units.Calorie;
import etomica.units.Mole;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.UnitRatio;

public class MethanolPotentialHelper {
	
	
	public static void initPotential(Space space, SpeciesMethanol species, PotentialGroup U_a_b, boolean pointCharges, double sigmaOC, double sigmaOH) {
		
		IPotentialAtomic u_O_O;
		IPotentialAtomic u_O_aC;
		IPotentialAtomic u_O_aH;
		IPotentialAtomic u_O_H;
	        
		IPotentialAtomic u_aC_aC;
		IPotentialAtomic u_aC_aH;
		IPotentialAtomic u_aC_H;
	        
		IPotentialAtomic u_aH_aH;
		IPotentialAtomic u_aH_H;
		IPotentialAtomic u_aH_X;
	        
		IPotentialAtomic u_H_H;

		P2RepRowley u_X_X;
	    
		/* 
		****************************************************************************
		****************************************************************************
	    Parameters for site-site interactions on different methanol molecules
	          
	    The original units are:
	          
	        kcal/mol for epsilon(i,j) - not in simulation units
	       	1/Angstroms for A(i,j) - in simulation units
	       	Angstroms for r and re - in simulation units
	        	  	
	    ****************************************************************************
        ****************************************************************************
        */
		
		// Create conversion factor to change epsilon values from kcal/mol to simulation units
		UnitRatio eunit = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
		
		if(pointCharges) {
			
			// oxygen and oxygen (O-O)
	        double epsilon_O_O = eunit.toSim(0.28128);
	        double A_O_O = 2.01169;
	        double re_O_O = 3.14862;
	        
	        // oxygen and the "alpha" carbon (O-aC)
	        double epsilon_O_aC = eunit.toSim(0.00049);
	        double A_O_aC = 2.87970;
	        double re_O_aC = 3.98955;
	        
	        // oxygen and the "alpha" hydrogen (O-aH)
	        double epsilon_O_aH = eunit.toSim(0.00953);
	        double A_O_aH = 2.39278;
	        double re_O_aH = 0.96556;
	        
	        // oxygen and hydrogen (O-H)
	        double epsilon_O_H = eunit.toSim(0.00000479);
	        double A_O_H = 1.44640;
	        double re_O_H = 6.75499;
	        
	        // "alpha" carbon and "alpha" carbon(aC-aC)
	        double epsilon_aC_aC = eunit.toSim(0.53431);
	        double A_aC_aC = 5.42872;
	        double re_aC_aC = 0.54673;
	        
	        // "alpha" carbon and "alpha" hydrogen (aC-aH)
	        double epsilon_aC_aH = eunit.toSim(1.20290);
	        double A_aC_aH = 1.48559;
	        double re_aC_aH = 2.41005;
	        
	        // "alpha" carbon and hydrogen (aC-H)
	        double epsilon_aC_H = eunit.toSim(0.30649);
	        double A_aC_H = 1.80093;
	        double re_aC_H = 2.67419;
	        
	        // "alpha" hydrogen and "alpha" hydrogen (aH-aH)
	        double epsilon_aH_aH = eunit.toSim(0.00303);
	        double A_aH_aH = 1.59429;
	        double re_aH_aH = 3.81180;
	        
	        // "alpha" hydrogen and hydrogen (aH-H)
	        double epsilon_aH_H = eunit.toSim(0.0000908);
	        double A_aH_H = 0.18342;
	        double re_aH_H = 1.38361;
	        
	        // fudge site and fudge site (X-X)
	        // Potential simplified to purely repulsive form: uXX = BXX * exp(-CXX*rXX)
	        double BXX = eunit.toSim(0.28457);
	        double CXX = 6.79460;
	        double rOX = 0.09310; // just used to draw X site
	        
	        // "alpha" hydrogen and fudge site (aH-X)
	        double epsilon_aH_X = eunit.toSim(0.03155);
	        double A_aH_X = 1.81525;
	        double re_aH_X = 2.98180;
	        
	        // hydrogen and hydrogen (H-H)
	        double epsilon_H_H = eunit.toSim(0.01048);
	        double A_H_H = 1.26072;
	        double re_H_H = 3.97536;
	        
	        // Point charges
	        double z_O = -0.6423;
	        double z_aC = 0.2550; // adapted from published value of 0.0001 s.t. molecule is neutral
	        double z_aH = 0.3873;
	        
	        /* 
	         * The point charges for the non-alpha hydrogens and the satellite site are zero: 
	         * One can use the unmodified Morse potential for interactions involving these sites.
	         */
	        
	        /* Cut-off distances for Coulombic interaction energies
	         * 
	         *    Because the term has the separation distance in the denominator, the coulombic interaction
	         *    can be enormously (and aphysically) attractive at very small (e.g. overlap) distances for sites having charges
	         *    of the same sign.
	         *    
	         *    For sites with charges of the same sign, a minimum distance is not required - the coulombic interaction 
	         *    effectively becomes a hard-core potential at very small distances.
	         */
	        
	        
	        /*
	        ****************************************************************************
	        ****************************************************************************
	        Directives for calculation of site-site interaction energies:
	        ****************************************************************************
	        ****************************************************************************
	        */
	        
	        u_O_O   = new P2ModifiedMorse(space, epsilon_O_O,   re_O_O,   A_O_O,   z_O,  z_O , 0   );
	        u_O_aC  = new P2ModifiedMorse(space, epsilon_O_aC,  re_O_aC,  A_O_aC,  z_O,  z_aC, sigmaOC  );
	        u_O_aH  = new P2ModifiedMorse(space, epsilon_O_aH,  re_O_aH,  A_O_aH,  z_O,  z_aH, sigmaOH  );
	        u_O_H   = new         P2Morse(space, epsilon_O_H,   re_O_H,   A_O_H   );
	        
	        u_aC_aC = new P2ModifiedMorse(space, epsilon_aC_aC, re_aC_aC, A_aC_aC, z_aC, z_aC, 0 );
	        u_aC_aH = new P2ModifiedMorse(space, epsilon_aC_aH, re_aC_aH, A_aC_aH, z_aC, z_aH, 0 );
	        u_aC_H  = new         P2Morse(space, epsilon_aC_H,  re_aC_H,  A_aC_H  );
	        
	        u_aH_aH = new P2ModifiedMorse(space, epsilon_aH_aH, re_aH_aH, A_aH_aH, z_aH, z_aH, 0 );
	        u_aH_H  = new         P2Morse(space, epsilon_aH_H,  re_aH_H,  A_aH_H  );
	        u_aH_X  = new         P2Morse(space, epsilon_aH_X,  re_aH_X,  A_aH_X  );
	        
	        u_H_H   = new         P2Morse(space, epsilon_H_H,   re_H_H,   A_H_H   );

	        u_X_X = new P2RepRowley (space, BXX, CXX);

		}
	        
		else {
			// oxygen and oxygen (O-O)
			double epsilon_O_O = eunit.toSim(0.56237);
			double A_O_O = 1.18149;
			double re_O_O = 3.74311;
			
			// oxygen and the "alpha" carbon (O-aC)
			double epsilon_O_aC = eunit.toSim(0.14519);
			double A_O_aC = 0.93092;
			double re_O_aC = 4.16094;
			
			// oxygen and the "alpha" hydrogen (O-aH)
			double epsilon_O_aH = eunit.toSim(12.46392);
			double A_O_aH = 1.62717;
			double re_O_aH = 1.50634;
			
			// oxygen and hydrogen (O-H)
			double epsilon_O_H = eunit.toSim(0.52961);
			double A_O_H = 1.28496;
			double re_O_H = 2.94761;
			
			// "alpha" carbon and "alpha" carbon(aC-aC)
			double epsilon_aC_aC = eunit.toSim(0.27110);
			double A_aC_aC = 3.18382;
			double re_aC_aC = 3.21345;
			
			// "alpha" carbon and "alpha" hydrogen (aC-aH)
			double epsilon_aC_aH = eunit.toSim(6.49019);
			double A_aC_aH = 12.24823;
			double re_aC_aH = 0.34630;
			
			// "alpha" carbon and hydrogen (aC-H)
			double epsilon_aC_H = eunit.toSim(0.42767);
			double A_aC_H = 4.77747;
			double re_aC_H = 0.50033;
		        
			// "alpha" hydrogen and "alpha" hydrogen (aH-aH)
			double epsilon_aH_aH = eunit.toSim(0.0000134);
			double A_aH_aH = 0.71522;
			double re_aH_aH = 11.27744;
			
			// "alpha" hydrogen and hydrogen (aH-H)
			double epsilon_aH_H = eunit.toSim(0.0000678);
			double A_aH_H = 0.38079;
			double re_aH_H = 14.70203;
		        
			// fudge site and fudge site (X-X)
			// Potential simplified to purely repulsive form: uXX = BXX * exp(-CXX*rXX)
			double BXX = eunit.toSim(14.63380);
			double CXX = 0.71468;
			double rOX = 0.96092; // just used to draw X site
		        
			// "alpha" hydrogen and fudge site (aH-X)// The satellite site, X, is closer to the oxygen atom in the model with point charges.
			double epsilon_aH_X = eunit.toSim(0.86651);
			double A_aH_X = 0.53543;
			double re_aH_X = 1.63407;
		        
			// hydrogen and hydrogen (H-H)
			double epsilon_H_H = eunit.toSim(0.01048);
			double A_H_H = 1.26072;
			double re_H_H = 3.97536;
			
			/*
			****************************************************************************
			****************************************************************************
		    Directives for calculation of site-site interaction energies:
			****************************************************************************
			****************************************************************************
			*/
		        
			u_O_O   = new P2Morse(space, epsilon_O_O,   re_O_O,   A_O_O   );
			u_O_aC  = new P2Morse(space, epsilon_O_aC,  re_O_aC,  A_O_aC  );
			u_O_aH  = new P2Morse(space, epsilon_O_aH,  re_O_aH,  A_O_aH  );
			u_O_H   = new P2Morse(space, epsilon_O_H,   re_O_H,   A_O_H   );
			
			u_aC_aC = new P2Morse(space, epsilon_aC_aC, re_aC_aC, A_aC_aC );
			u_aC_aH = new P2Morse(space, epsilon_aC_aH, re_aC_aH, A_aC_aH );
			u_aC_H  = new P2Morse(space, epsilon_aC_H,  re_aC_H,  A_aC_H  );
			
			u_aH_aH = new P2Morse(space, epsilon_aH_aH, re_aH_aH, A_aH_aH );
			u_aH_H  = new P2Morse(space, epsilon_aH_H,  re_aH_H,  A_aH_H  );
			u_aH_X  = new P2Morse(space, epsilon_aH_X,  re_aH_X,  A_aH_X  );
		        
			u_H_H   = new P2Morse(space, epsilon_H_H,   re_H_H,   A_H_H   );
	
			u_X_X = new P2RepRowley (space, BXX, CXX);
		}
		
		
	        
		IAtomType type_O  = species.getOxygenType();
		IAtomType type_aC = species.getAlphaCarbonType(); 
		IAtomType type_aH = species.getAlphaHydrogenType();
		IAtomType type_H  = species.getHydrogenType();
		IAtomType type_X  = species.getXType();
	        
		/*
		****************************************************************************
		****************************************************************************
	    Directives for calculation of the MOLECULAR potential, U_a_b
		****************************************************************************
		****************************************************************************
		*/

	         
		U_a_b.addPotential(u_O_O,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_O }));
		
		U_a_b.addPotential(u_O_aC,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_aC}));
		U_a_b.addPotential(u_O_aC,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC, type_O }));
	         
		U_a_b.addPotential(u_O_aH,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_aH}));
		U_a_b.addPotential(u_O_aH,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH, type_O }));
	         
		U_a_b.addPotential(u_O_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_O,  type_H }));
		U_a_b.addPotential(u_O_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,  type_O }));
	         
	         
		U_a_b.addPotential(u_aC_aC,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_aC}));
	         
		U_a_b.addPotential(u_aC_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_aH}));
		U_a_b.addPotential(u_aC_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_aC}));
	         
		U_a_b.addPotential(u_aC_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aC,  type_H }));
		U_a_b.addPotential(u_aC_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_aC}));
	         
	         
		U_a_b.addPotential(u_aH_aH,  ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_aH}));
	         
		U_a_b.addPotential(u_aH_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_H }));
		U_a_b.addPotential(u_aH_H,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_aH}));
	         
		U_a_b.addPotential(u_aH_X,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_aH,  type_X }));
		U_a_b.addPotential(u_aH_X,   ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_X,   type_aH}));
	         
	         
		U_a_b.addPotential(u_H_H,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_H,   type_H }));
	         
	         
		U_a_b.addPotential(u_X_X,    ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{type_X,   type_X }));
		
	}
	

}
