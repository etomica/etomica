/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.api.IPotentialAtomic;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.potential.P2ModifiedMorse;
import etomica.potential.P2Morse;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.units.*;

public class EthanolPotentialHelper {
	
	public static void initPotential(Space space, SpeciesEthanol species, PotentialGroup U_a_b, boolean pointCharges) {
		
		IPotentialAtomic u_O_O;
		IPotentialAtomic u_O_aC;
		IPotentialAtomic u_O_C;
		IPotentialAtomic u_O_aH;
		IPotentialAtomic u_O_H;
	        
		IPotentialAtomic u_aC_aC;
		IPotentialAtomic u_aC_C;
		IPotentialAtomic u_aC_aH;
		IPotentialAtomic u_aC_H;
		
		IPotentialAtomic u_C_aH;
	        
		IPotentialAtomic u_aH_aH;
		IPotentialAtomic u_aH_H;
		IPotentialAtomic u_aH_X;
		
		IPotentialAtomic u_C_C;
		IPotentialAtomic u_C_H;
	        
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
	        double epsilon_O_O = eunit.toSim(8.49412);
	        double a_O_O = 0.92972;
	        double re_O_O = 0.67646;
	        
	        // oxygen and the "alpha" carbon (O-aC)
	        double epsilon_O_aC = eunit.toSim(0.0000179);
	        double a_O_aC = 0.77318;
	        double re_O_aC = 2.60416;
	        
	        // oxygen and the non-alpha carbon (O-C)
	        double epsilon_O_C = eunit.toSim(7.16564);
	        double a_O_C = 1.46978;
	        double re_O_C = 2.57620;
	        
	        // oxygen and the "alpha" hydrogen (O-aH)
	        double epsilon_O_aH = eunit.toSim(3.61543);
	        double a_O_aH = 1.85755;
	        double re_O_aH = 1.93924;
	        
	        // oxygen and hydrogen (O-H)
	        double epsilon_O_H = eunit.toSim(0.0000617);
	        double a_O_H = 0.94407;
	        double re_O_H = 8.28351;
	        
	        // "alpha" carbon and "alpha" carbon(aC-aC)
	        double epsilon_aC_aC = eunit.toSim(4.59505);
	        double a_aC_aC = 1.29692;
	        double re_aC_aC = 2.92735;
	        
	        // "alpha" carbon and non-alpha carbon(aC-C)
	        double epsilon_aC_C = eunit.toSim(0.00978);
	        double a_aC_C = 0.32849;
	        double re_aC_C = 10.9008;
	        
	        // "alpha" carbon and "alpha" hydrogen (aC-aH)
	        double epsilon_aC_aH = eunit.toSim(1.74991);
	        double a_aC_aH = 2.82282;
	        double re_aC_aH = 2.38571;
	        
	        // "alpha" carbon and hydrogen (aC-H)
	        double epsilon_aC_H = eunit.toSim(5.72187);
	        double a_aC_H = 4.84533;
	        double re_aC_H = 0.49777;
	        
	        // non-alpha carbon and "alpha" hydrogen (C-aH)
	        double epsilon_C_aH = eunit.toSim(0.06512);
	        double a_C_aH = 10.01582;
	        double re_C_aH = 0.05860;
	        
	        // "alpha" hydrogen and "alpha" hydrogen (aH-aH)
	        double epsilon_aH_aH = eunit.toSim(0.01461);
	        double a_aH_aH = 0.94465;
	        double re_aH_aH = 5.00333;
	        
	        // "alpha" hydrogen and hydrogen (aH-H)
	        double epsilon_aH_H = eunit.toSim(0.68571);
	        double a_aH_H = 4.62914;
	        double re_aH_H = 0.50073;
	        
	        // fudge site and fudge site (X-X)
	        // Potential simplified to purely repulsive form: uXX = BXX * exp(-CXX*rXX)
	        double BXX = eunit.toSim(1.00665);
	        double CXX = 4.24048;
	        double rOX = 0.33858; // just used to draw X site
	        
	        // "alpha" hydrogen and fudge site (aH-X)// The satellite site, X, is closer to the oxygen atom in the model with point charges.
	        double epsilon_aH_X = eunit.toSim(0.00105);
	        double a_aH_X = 0.66751;
	        double re_aH_X = 7.88000;
	        
	        // non-alpha carbon and non-alpha carbon (C-C)
	        double epsilon_C_C = eunit.toSim(0.05133);
	        double a_C_C = 1.45985;
	        double re_C_C = 4.34117;
	        
	        // non-alpha carbon and hydrogen (C-H)
	        double epsilon_C_H = eunit.toSim(0.35562);
	        double a_C_H = 2.11174;
	        double re_C_H = 2.60211;
	        
	        // hydrogen and hydrogen (H-H)
	        double epsilon_H_H = eunit.toSim(0.01048);
	        double a_H_H = 1.26072;
	        double re_H_H = 3.97536;
	        
	        // Point charges
	        double z_O = -0.7217;
	        double z_aC = 0.3318;
	        double z_aH = 0.3899;
	        
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
	        
	        double cc_O_O  = 0;
	        double cc_O_aC = 0.5;
	        double cc_O_aH = 0.3;
	        	
	    	double cc_aC_aC = 0;
	    	double cc_aC_aH = 0;
	        	
	    	double cc_aH_aH = 0;
	        
	        /*
	        ****************************************************************************
	        ****************************************************************************
	        Directives for calculation of site-site interaction energies:
	        ****************************************************************************
	        ****************************************************************************
	        */
	        
	        u_O_O   = new P2ModifiedMorse(space, epsilon_O_O,   re_O_O,   a_O_O,   z_O,  z_O , cc_O_O   );
	        u_O_aC  = new P2ModifiedMorse(space, epsilon_O_aC,  re_O_aC,  a_O_aC,  z_O,  z_aC, cc_O_aC  );
	        u_O_C   = new P2Morse        (space, epsilon_O_C,   re_O_C,   a_O_C   );
	        u_O_aH  = new P2ModifiedMorse(space, epsilon_O_aH,  re_O_aH,  a_O_aH,  z_O,  z_aH, cc_O_aH  );
	        u_O_H   = new P2Morse        (space, epsilon_O_H,   re_O_H,   a_O_H   );
	        
	        u_aC_aC = new P2ModifiedMorse(space, epsilon_aC_aC, re_aC_aC, a_aC_aC, z_aC, z_aC, cc_aC_aC );
	        u_aC_C  = new P2Morse        (space, epsilon_aC_C,  re_aC_C,  a_aC_C  );
	        u_aC_aH = new P2ModifiedMorse(space, epsilon_aC_aH, re_aC_aH, a_aC_aH, z_aC, z_aH, cc_aC_aH );
	        u_aC_H  = new P2Morse        (space, epsilon_aC_H,  re_aC_H,  a_aC_H  );
	        
	        u_C_aH  = new P2Morse        (space, epsilon_C_aH,  re_C_aH,  a_C_aH  );
	        
	        u_aH_aH = new P2ModifiedMorse(space, epsilon_aH_aH, re_aH_aH, a_aH_aH, z_aH, z_aH, cc_aH_aH );
	        u_aH_H  = new P2Morse        (space, epsilon_aH_H,  re_aH_H,  a_aH_H  );
	        u_aH_X  = new P2Morse        (space, epsilon_aH_X,  re_aH_X,  a_aH_X  );
	        
	        u_C_C   = new P2Morse        (space, epsilon_C_C,   re_C_C,   a_C_C   );
	        u_C_H   = new P2Morse        (space, epsilon_C_H,   re_C_H,   a_C_H   );
	        
	        u_H_H   = new P2Morse        (space, epsilon_H_H,   re_H_H,   a_H_H   );

	        u_X_X = new P2RepRowley (space, BXX, CXX);

		}
		else {
			
			// oxygen and oxygen (O-O)
	        double epsilon_O_O = eunit.toSim(0.03292);
	        double a_O_O = 0.89878;
	        double re_O_O = 5.8420;
	        
	        // oxygen and the "alpha" carbon (O-aC)
	        double epsilon_O_aC = eunit.toSim(2.31319);
	        double a_O_aC = 1.84553;
	        double re_O_aC = 2.99986;
	        
	        // oxygen and the non-alpha carbon (O-C)
	        double epsilon_O_C = eunit.toSim(2.62034);
	        double a_O_C = 1.79601;
	        double re_O_C = 2.95039;
	        
	        // oxygen and the "alpha" hydrogen (O-aH)
	        double epsilon_O_aH = eunit.toSim(18.0197);
	        double a_O_aH = 1.64599;
	        double re_O_aH = 1.47544;
	        
	        // oxygen and hydrogen (O-H)
	        double epsilon_O_H = eunit.toSim(0.00254);
	        double a_O_H = 1.01726;
	        double re_O_H = 5.77252;
	        
	        // "alpha" carbon and "alpha" carbon(aC-aC)
	        double epsilon_aC_aC = eunit.toSim(0.00077);
	        double a_aC_aC = 6.09818;
	        double re_aC_aC = 3.78185;
	        
	        // "alpha" carbon and non-alpha carbon(aC-C)
	        double epsilon_aC_C = eunit.toSim(0.18367);
	        double a_aC_C = 1.66760;
	        double re_aC_C = 3.82351;
	        
	        // "alpha" carbon and "alpha" hydrogen (aC-aH)
	        double epsilon_aC_aH = eunit.toSim(6.34217);
	        double a_aC_aH = 12.6701;
	        double re_aC_aH = 0.23124;
	        
	        // "alpha" carbon and hydrogen (aC-H)
	        double epsilon_aC_H = eunit.toSim(6.92482);
	        double a_aC_H = 9.17206;
	        double re_aC_H = 0.16668;
	        
	        // non-alpha carbon and "alpha" hydrogen (C-aH)
	        double epsilon_C_aH = eunit.toSim(0.07353);
	        double a_C_aH = 5.32819;
	        double re_C_aH = 0.00775;
	        
	        // "alpha" hydrogen and "alpha" hydrogen (aH-aH)
	        double epsilon_aH_aH = eunit.toSim(0.0000628);
	        double a_aH_aH = 0.70122;
	        double re_aH_aH = 10.48771;
	        
	        // "alpha" hydrogen and hydrogen (aH-H)
	        double epsilon_aH_H = eunit.toSim(0.000000962);
	        double a_aH_H = 1.50906;
	        double re_aH_H = 6.41108;
	        
	        // fudge site and fudge site (X-X)
	        // Potential simplified to purely repulsive form: uXX = BXX * exp(-CXX*rXX)
	        double BXX = eunit.toSim(5.99551);
	        double CXX = 0.60331;
	        double rOX = 1.10795; // just used to draw X site
	        
	        // "alpha" hydrogen and fudge site (aH-X)// The satellite site, X, is closer to the oxygen atom in the model with point charges.
	        double epsilon_aH_X = eunit.toSim(0.70186);
	        double a_aH_X = 0.84393;
	        double re_aH_X = 1.90702;
	        
	        // non-alpha carbon and non-alpha carbon (C-C)
	        double epsilon_C_C = eunit.toSim(0.05133);
	        double a_C_C = 1.45985;
	        double re_C_C = 4.34117;
	        
	        // non-alpha carbon and hydrogen (C-H)
	        double epsilon_C_H = eunit.toSim(0.35562);
	        double a_C_H = 2.11174;
	        double re_C_H = 2.60211;
	        
	        // hydrogen and hydrogen (H-H)
	        double epsilon_H_H = eunit.toSim(0.01048);
	        double a_H_H = 1.26072;
	        double re_H_H = 3.97536;
	        
	        /*
	        ****************************************************************************
	        ****************************************************************************
	        Directives for calculation of site-site interaction energies:
	        ****************************************************************************
	        ****************************************************************************
	        */
	        
	        u_O_O   = new P2Morse(space, epsilon_O_O,   re_O_O,   a_O_O   );
	        u_O_aC  = new P2Morse(space, epsilon_O_aC,  re_O_aC,  a_O_aC  );
	        u_O_C   = new P2Morse(space, epsilon_O_C,   re_O_C,   a_O_C   );
	        u_O_aH  = new P2Morse(space, epsilon_O_aH,  re_O_aH,  a_O_aH  );
	        u_O_H   = new P2Morse(space, epsilon_O_H,   re_O_H,   a_O_H   );
	        
	        u_aC_aC = new P2Morse(space, epsilon_aC_aC, re_aC_aC, a_aC_aC );
	        u_aC_C  = new P2Morse(space, epsilon_aC_C,  re_aC_C,  a_aC_C  );
	        u_aC_aH = new P2Morse(space, epsilon_aC_aH, re_aC_aH, a_aC_aH );
	        u_aC_H  = new P2Morse(space, epsilon_aC_H,  re_aC_H,  a_aC_H  );
	        
	        u_C_aH  = new P2Morse(space, epsilon_C_aH,  re_C_aH,  a_C_aH  );
	        
	        u_aH_aH = new P2Morse(space, epsilon_aH_aH, re_aH_aH, a_aH_aH );
	        u_aH_H  = new P2Morse(space, epsilon_aH_H,  re_aH_H,  a_aH_H  );
	        u_aH_X  = new P2Morse(space, epsilon_aH_X,  re_aH_X,  a_aH_X  );
	        
	        u_C_C   = new P2Morse(space, epsilon_C_C,   re_C_C,   a_C_C   );
	        u_C_H   = new P2Morse(space, epsilon_C_H,   re_C_H,   a_C_H   );
	        
	        u_H_H   = new P2Morse(space, epsilon_H_H,   re_H_H,   a_H_H   );

	        u_X_X = new P2RepRowley (space, BXX, CXX);

	        
		}
		
		/*
		****************************************************************************
		****************************************************************************
	    Create instances of the types of molecular sites
		****************************************************************************
		****************************************************************************
		*/


        AtomType type_O = species.getOxygenType();
        AtomType type_aC = species.getAlphaCarbonType();
        AtomType type_C = species.getCarbonType();
        AtomType type_aH = species.getAlphaHydrogenType();
        AtomType type_H = species.getHydrogenType();
        AtomType type_X = species.getXType();

		/*
		****************************************************************************
		****************************************************************************
	    Directives for calculation of the MOLECULAR potential, U_a_b
		****************************************************************************
		****************************************************************************
		*/


        U_a_b.addPotential(u_O_O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_O, type_O}));

        U_a_b.addPotential(u_O_aC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_O, type_aC}));
        U_a_b.addPotential(u_O_aC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aC, type_O}));

        U_a_b.addPotential(u_O_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_O, type_C}));
        U_a_b.addPotential(u_O_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_C, type_O}));

        U_a_b.addPotential(u_O_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_O, type_aH}));
        U_a_b.addPotential(u_O_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_O}));

        U_a_b.addPotential(u_O_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_O, type_H}));
        U_a_b.addPotential(u_O_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_H, type_O}));


        U_a_b.addPotential(u_aC_aC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aC, type_aC}));

        U_a_b.addPotential(u_aC_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aC, type_C}));
        U_a_b.addPotential(u_aC_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_C, type_aC}));

        U_a_b.addPotential(u_aC_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aC, type_aH}));
        U_a_b.addPotential(u_aC_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_aC}));

        U_a_b.addPotential(u_aC_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aC, type_H}));
        U_a_b.addPotential(u_aC_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_H, type_aC}));


        U_a_b.addPotential(u_C_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_C, type_aH}));
        U_a_b.addPotential(u_C_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_C}));


        U_a_b.addPotential(u_aH_aH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_aH}));

        U_a_b.addPotential(u_aH_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_H}));
        U_a_b.addPotential(u_aH_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_H, type_aH}));

        U_a_b.addPotential(u_aH_X, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_aH, type_X}));
        U_a_b.addPotential(u_aH_X, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_X, type_aH}));


        U_a_b.addPotential(u_C_C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_C, type_C}));

        U_a_b.addPotential(u_C_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_C, type_H}));
        U_a_b.addPotential(u_C_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_H, type_C}));


        U_a_b.addPotential(u_H_H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_H, type_H}));


        U_a_b.addPotential(u_X_X, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{type_X, type_X}));


    }

}
