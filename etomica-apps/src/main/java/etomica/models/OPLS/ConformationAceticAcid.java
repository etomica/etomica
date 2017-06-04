/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for acetic acid
 * 
 * @author Hye Min Kim
 * Nov, 2011
 */
public class ConformationAceticAcid implements IConformation {

    public ConformationAceticAcid(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
    	
    	double bondCH3C = 1.52;//Angstrom, (CH3)-C
    	double bondCSBO = 1.364;//Angstrom, C-O
    	double bondCDBO = 1.214;//Angstrom, C=O
    	double bondOH = 0.97;//Angstrom, O-H
    	double bondCH = 1.890782497;//Angstrom, C-H
    	double angleCCDBO = 126*Math.PI/180;//C-C=O
    	double angleOCO = 123*Math.PI/180;//O=C-O
    	double angleCOH = 107*Math.PI/180;//C-O-H
    	double angleCCSBO =111*Math.PI/180;//C-C-O
    	double angleCCH =140.37995*Math.PI/180;//C-C-H
    	
    	IAtom cH3 = list.getAtom(SpeciesAceticAcid.indexCH3);
        cH3.getPosition().E(new double[] {bondCH3C, 0.0, 0.0});
        
        IAtom c = list.getAtom(SpeciesAceticAcid.indexC);
        c.getPosition().E(new double[] {0.0, 0.0, 0.0});
        
        IAtom dBO = list.getAtom(SpeciesAceticAcid.indexDBO);
        dBO.getPosition().E(new double[] {bondCDBO*Math.cos(angleCCDBO), bondCDBO*Math.sin(angleCCDBO), 0.0});
        
        IAtom sBO = list.getAtom(SpeciesAceticAcid.indexSBO);
        sBO.getPosition().E(new double[] {bondCSBO*Math.cos(-angleCCSBO), bondCSBO*Math.sin(-angleCCSBO), 0.0});//opposite direction
        
        IAtom hydrogen = list.getAtom(SpeciesAceticAcid.indexH);
        hydrogen.getPosition().E(new double[] {bondCH*Math.cos(-angleCCH), bondCH*Math.sin(-angleCCH), 0.0});//opposite direction       
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;

}

