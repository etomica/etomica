/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.traPPE;

import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.potential.P2ElectrostaticWithHardCore;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Constants;

/*
 * PotentialHelper class for TraPPE model of methanol: http://www.chem.umn.edu/groups/siepmann/trappe/intro.php
 * 
 * CH3 group of methanol is united into a single site 
 * 
 * K.R. Schadel 2008
 */

public class MethanolPotentialHelper {
	
	
	public static void initPotential(Space space, SpeciesGeneral species, PotentialGroup U_a_b) {
		

	        
		
		/* 
		****************************************************************************
		****************************************************************************
	    Parameters for site-site interactions on different methanol molecules
	          
	    The original units are:
	          
	        K for epsilon - really is epsilon/kB
	       	Angstroms for sigma
	        	  	
	    ****************************************************************************
        ****************************************************************************
        */
		
		// TraPPE params given in epsilon/kB [=] K;
		// Remember that Constants.BOLTZMANN_K is in simulation units
		double epsilonCH3 = Kelvin.UNIT.toSim(98.00); 
		double epsilonO = Kelvin.UNIT.toSim(93.00); 
		//double epsilonH = 0; // There is  no LJ interaction involving the hydrogens!
		
		// sigma [=] Angstroms;
		double sigmaCH3 = 3.75;
		double sigmaO = 3.02;
		double sigmaH = 0;
		
		//Lorentz-Berthelot combining rules:
		
		double epsilonCH3O = Math.sqrt(epsilonCH3*epsilonO);
		//double epsilonCH3H = Math.sqrt(epsilonCH3*epsilonH);
		//double epsilonOH   = Math.sqrt(epsilonO*epsilonH);
		
		/*epsilonCH3 = Kelvin.UNIT.toSim(epsilonCH3);
		epsilonO = Kelvin.UNIT.toSim(epsilonO);
		epsilonH = Kelvin.UNIT.toSim(epsilonH);
		epsilonCH3O = Kelvin.UNIT.toSim(epsilonCH3O);
		epsilonCH3H = Kelvin.UNIT.toSim(epsilonCH3H);
		epsilonOH = Kelvin.UNIT.toSim(epsilonOH);*/
		
		double sigmaCH3O = 0.5*(sigmaCH3+sigmaO);
		double sigmaCH3H = 0.5*(sigmaCH3+sigmaH);
		double sigmaOH   = 0.5*(sigmaO  +sigmaH);
			 
        // Point charges
		
		double hatedConversionTerm = Math.sqrt(4*Math.PI*Constants.EPSILON_0);
		
		/*double zCH3 =  0.265/hatedConversionTerm;
        double zO   = -0.700/hatedConversionTerm;
        double zH   =  0.435/hatedConversionTerm;*/
		
		double zCH3 =  Electron.UNIT.toSim(0.265);
        double zO   = Electron.UNIT.toSim(-0.7);
        double zH   =  Electron.UNIT.toSim(0.435);
        
        
	          
	        
        /*
        ****************************************************************************
        ****************************************************************************
        Directives for calculation of site-site interaction energies:
        ****************************************************************************
        ****************************************************************************
        */
        
        // LJ site-site interaction energies
        
        P2LennardJones uLJCH3CH3 = new P2LennardJones(space);
        P2LennardJones uLJCH3O   = new P2LennardJones(space);     
        P2LennardJones uLJOO     = new P2LennardJones(space);

        uLJCH3CH3.setSigma(sigmaCH3);
        uLJCH3O.setSigma(sigmaCH3O);
        uLJOO.setSigma(sigmaO);
        
        uLJCH3CH3.setEpsilon(epsilonCH3);
        uLJCH3O.setEpsilon(epsilonCH3O);
        uLJOO.setEpsilon(epsilonO);
        
        // Coulombic site-site interaction energies
        
        P2ElectrostaticWithHardCore uCH3CH3 = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3O   = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3H   = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uOO     = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uOH     = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uHH     = new P2ElectrostaticWithHardCore(space);
        
        uCH3CH3.setCharge1(zCH3);
        uCH3CH3.setCharge2(zCH3);
        uCH3CH3.setSigma(0);
        
        uCH3O.setCharge1(zCH3);
        uCH3O.setCharge2(zO);
        uCH3O.setSigma(0.1);
        // sigmaCH3O = 3.385 A
        // 0.7 * sigmaCH30 = 2.3695 A
        
        uCH3H.setCharge1(zCH3);
        uCH3H.setCharge2(zH);
        uCH3H.setSigma(0);
        
        uOO.setCharge1(zO);
        uOO.setCharge2(zO);
        uOO.setSigma(0);
        // sigmaO = 3.02
        // H can come within 0.7*3.02 - 0.95 = 1.164 A of O.  
        
        uOH.setCharge1(zO);
        uOH.setCharge2(zH);
        uOH.setSigma(0.1);
        // sigmaO = 3.02;
        // sigmaOH = 1.51;
        // 0.7 * sigmaO = 2.114
        // 0.7 * sigmaO + bondOH = 3.064
        
        uHH.setCharge1(zH);
        uHH.setCharge2(zH);
        uHH.setSigma(0);


        AtomType typeCH3 = species.getTypeByName("CH3");
        AtomType typeO = species.getTypeByName("O");
        AtomType typeH = species.getTypeByName("H");

		/*
		****************************************************************************
		****************************************************************************
	    Directives for calculation of the MOLECULAR potential, U_a_b
		****************************************************************************
		****************************************************************************
		*/

        U_a_b.addPotential(uLJCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));

        U_a_b.addPotential(uLJCH3O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeO}));
        U_a_b.addPotential(uLJCH3O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeCH3}));

        U_a_b.addPotential(uLJOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));


        U_a_b.addPotential(uCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));

        U_a_b.addPotential(uCH3O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeO}));
        U_a_b.addPotential(uCH3O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeCH3}));

        U_a_b.addPotential(uCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeH}));
        U_a_b.addPotential(uCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH3}));

        U_a_b.addPotential(uOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));

        U_a_b.addPotential(uOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeH}));
        U_a_b.addPotential(uOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeO}));


        U_a_b.addPotential(uHH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));

    }

}

