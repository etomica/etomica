/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.potential.P2ElectrostaticWithHardCore;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * PotentialHelper class of acetic acid using OPLS united atom model
 * 
 * @author Hye Min Kim
 * Dec, 2011
 */
public class AceticAcidPotentialHelper {
	
	
	public static void initPotential(Space space, SpeciesAceticAcid species, PotentialGroup p) {
		double epsilonCH3 = Kelvin.UNIT.toSim(80.6); 
		double epsilonC = Kelvin.UNIT.toSim(52.9);
		double epsilonDBO = Kelvin.UNIT.toSim(105.7);
		double epsilonSBO = Kelvin.UNIT.toSim(85.6);
		double epsilonH = 0;
		double sigmaCH3 = 3.91;
		double sigmaC = 3.75;
		double sigmaDBO = 2.96;
		double sigmaSBO = 3.00;
		double sigmaH = 0;
		
		double epsilonCH3C = Math.sqrt(epsilonCH3*epsilonC);//L-B combining rule
		double epsilonCH3DBO = Math.sqrt(epsilonCH3*epsilonDBO);
		double epsilonCH3SBO = Math.sqrt(epsilonCH3*epsilonSBO);
		double epsilonCDBO = Math.sqrt(epsilonC*epsilonDBO);
		double epsilonCSBO = Math.sqrt(epsilonC*epsilonSBO);
		double epsilonDBOSBO = Math.sqrt(epsilonDBO*epsilonSBO);
		
		double sigmaCH3C = 0.5*(sigmaCH3+sigmaC);
		double sigmaCH3DBO = 0.5*(sigmaCH3+sigmaDBO);
		double sigmaCH3SBO = 0.5*(sigmaCH3+sigmaSBO);
		double sigmaCDBO = 0.5*(sigmaC+sigmaDBO);
		double sigmaCSBO = 0.5*(sigmaC+sigmaSBO);
		double sigmaDBOSBO = 0.5*(sigmaDBO+sigmaSBO);
		
		double zCH3 =  Electron.UNIT.toSim(0.08);//partial charge of CH3 site
		double zC =  Electron.UNIT.toSim(0.55);
        double zDBO   = Electron.UNIT.toSim(-0.50);
        double zSBO   = Electron.UNIT.toSim(-0.58);
        double zH   =  Electron.UNIT.toSim(0.45);

        P2LennardJones uLJCH3CH3 = new P2LennardJones(space);
        P2LennardJones uLJCH3C = new P2LennardJones(space);
        P2LennardJones uLJCH3DBO = new P2LennardJones(space);
        P2LennardJones uLJCH3SBO = new P2LennardJones(space);  
        P2LennardJones uLJCC = new P2LennardJones(space);  
        P2LennardJones uLJCDBO = new P2LennardJones(space);
        P2LennardJones uLJCSBO = new P2LennardJones(space);  
        P2LennardJones uLJDBODBO = new P2LennardJones(space);
        P2LennardJones uLJDBOSBO = new P2LennardJones(space);
        P2LennardJones uLJSBOSBO = new P2LennardJones(space);

        uLJCH3CH3.setSigma(sigmaCH3);
        uLJCH3C.setSigma(sigmaCH3C);
        uLJCH3DBO.setSigma(sigmaCH3DBO);
        uLJCH3SBO.setSigma(sigmaCH3SBO);
        uLJCC.setSigma(sigmaC);
        uLJCDBO.setSigma(sigmaCDBO);
        uLJCSBO.setSigma(sigmaCSBO);
        uLJDBODBO.setSigma(sigmaDBO);
        uLJDBOSBO.setSigma(sigmaDBOSBO);
        uLJSBOSBO.setSigma(sigmaSBO);
        
        uLJCH3CH3.setEpsilon(epsilonCH3);
        uLJCH3C.setEpsilon(epsilonCH3C);
        uLJCH3DBO.setEpsilon(epsilonCH3DBO);
        uLJCH3SBO.setEpsilon(epsilonCH3SBO);
        uLJCC.setEpsilon(epsilonC);
        uLJCDBO.setEpsilon(epsilonCDBO);
        uLJCSBO.setEpsilon(epsilonCSBO);
        uLJDBODBO.setEpsilon(epsilonDBO);
        uLJDBOSBO.setEpsilon(epsilonDBOSBO);
        uLJSBOSBO.setEpsilon(epsilonSBO);
        
        // Coulombic site-site interaction energies
        
        P2ElectrostaticWithHardCore uCH3CH3 = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3C = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3DBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3SBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH3H = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uCC = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCDBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCSBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uCH = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uDBODBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uDBOSBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uDBOH = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uSBOSBO = new P2ElectrostaticWithHardCore(space);
        P2ElectrostaticWithHardCore uSBOH = new P2ElectrostaticWithHardCore(space);
        
        P2ElectrostaticWithHardCore uHH = new P2ElectrostaticWithHardCore(space);
        
        uCH3CH3.setCharge1(zCH3);
        uCH3CH3.setCharge2(zCH3);
        uCH3CH3.setSigma(0);
        
        uCH3C.setCharge1(zCH3);
        uCH3C.setCharge2(zC);
        uCH3C.setSigma(0);
        
        uCH3DBO.setCharge1(zCH3);
        uCH3DBO.setCharge2(zDBO);
        uCH3DBO.setSigma(0.1);
        
        uCH3SBO.setCharge1(zCH3);
        uCH3SBO.setCharge2(zSBO);
        uCH3SBO.setSigma(0.1);
        
        uCH3H.setCharge1(zCH3);
        uCH3H.setCharge2(zH);
        uCH3H.setSigma(0);
        
        uCC.setCharge1(zC);
        uCC.setCharge2(zC);
        uCC.setSigma(0);
        
        uCDBO.setCharge1(zC);
        uCDBO.setCharge2(zDBO);
        uCDBO.setSigma(0.1);
        
        uCSBO.setCharge1(zC);
        uCSBO.setCharge2(zSBO);
        uCSBO.setSigma(0.1);
        
        uCH.setCharge1(zC);
        uCH.setCharge2(zH);
        uCH.setSigma(0);
        
        uDBODBO.setCharge1(zDBO);
        uDBODBO.setCharge2(zDBO);
        uDBODBO.setSigma(0);
        
        uDBOSBO.setCharge1(zDBO);
        uDBOSBO.setCharge2(zSBO);
        uDBOSBO.setSigma(0);
        
        uDBOH.setCharge1(zDBO);
        uDBOH.setCharge2(zH);
        uDBOH.setSigma(0.1);
        
        uSBOSBO.setCharge1(zSBO);
        uSBOSBO.setCharge2(zSBO);
        uSBOSBO.setSigma(0);
        
        uSBOH.setCharge1(zSBO);
        uSBOH.setCharge2(zH);
        uSBOH.setSigma(0.1);

        uHH.setCharge1(zH);
        uHH.setCharge2(zH);
        uHH.setSigma(0);

        AtomType typeCH3 = species.getCH3Type();
        AtomType typeC = species.getCType();
        AtomType typeDBO = species.getDBOType();
        AtomType typeSBO = species.getSBOType();
        AtomType typeH = species.getHType();


        p.addPotential(uLJCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));

        p.addPotential(uLJCH3C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeC}));
        p.addPotential(uLJCH3C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH3}));

        p.addPotential(uLJCH3DBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeDBO}));
        p.addPotential(uLJCH3DBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeCH3}));

        p.addPotential(uLJCH3SBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeSBO}));
        p.addPotential(uLJCH3SBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeCH3}));

        p.addPotential(uLJCC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));

        p.addPotential(uLJCDBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeDBO}));
        p.addPotential(uLJCDBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeC}));

        p.addPotential(uLJCSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeSBO}));
        p.addPotential(uLJCSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeC}));

        p.addPotential(uLJDBODBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeDBO}));

        p.addPotential(uLJDBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeSBO}));
        p.addPotential(uLJDBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeDBO}));

        p.addPotential(uLJSBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeSBO}));

        p.addPotential(uCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));

        p.addPotential(uCH3C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeC}));
        p.addPotential(uCH3C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeCH3}));

        p.addPotential(uCH3DBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeDBO}));
        p.addPotential(uCH3DBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeCH3}));

        p.addPotential(uCH3SBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeSBO}));
        p.addPotential(uCH3SBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeCH3}));

        p.addPotential(uCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeH}));
        p.addPotential(uCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH3}));

        p.addPotential(uCC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));

        p.addPotential(uCDBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeDBO}));
        p.addPotential(uCDBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeC}));

        p.addPotential(uCSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeSBO}));
        p.addPotential(uCSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeC}));

        p.addPotential(uCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeH}));
        p.addPotential(uCH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeC}));

        p.addPotential(uDBODBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeDBO}));

        p.addPotential(uDBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeSBO}));
        p.addPotential(uDBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeDBO}));

        p.addPotential(uDBOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeDBO, typeH}));
        p.addPotential(uDBOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeDBO}));

        p.addPotential(uSBOSBO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeSBO}));

        p.addPotential(uSBOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeSBO, typeH}));
        p.addPotential(uSBOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeSBO}));

        p.addPotential(uHH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));

    }

}

