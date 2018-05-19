/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.space.Space;
import etomica.units.Degree;

/**
 * This reads from a configuration file containing only H, H and O.  The M
 * position is determined from the others.  Also, the geometry is relaxed due
 * to the file containing limited precision.
 */
public class ConfigurationFileTIP4P extends ConfigurationFile {
	protected final Space space;
	protected boolean isIce;
	public ConfigurationFileTIP4P(String aConfName, Space space, boolean isIce) {
		super(aConfName);
		this.space = space;
		this.isIce = isIce;
	}
    public void initializeCoordinates(Box box) {
        IAtomList leafList = box.getLeafList();
        String fileName = confName+".pos";
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            int nLeaf = leafList.size();
            Vector tmp = space.makeVector();
            Vector tmp2 = space.makeVector();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.get(iLeaf);
                if(a.getLeafIndex() % 4 != 3 ){// skip M (order is HHOM)
                	
                    //~ Ignore M positions
                	setPosition(a,bufReader.readLine());
                	//Added to solve the round-off error of having atoms at the edge
//                	a.getPosition().PE(0.001);
                }else{ // =3
                    //~ Define M positions
                	Vector h1 = leafList.get(iLeaf - 3).getPosition();//0 H1
                	Vector h2 = leafList.get(iLeaf - 2).getPosition();//1 H2
                	Vector o = leafList.get(iLeaf - 1).getPosition();//2 O
                	
					if(true){
	                	//Sice I/P xyz has ONLY 6 digits we do this reconstruction :)
	                	h1.ME(o);
	                	h1.normalize();
	                	h1.TE(0.9572);//Same for TIP4P & TIP4P/Ice
	                	h1.PE(o);
	                	tmp.Ev1Mv2(h1, o);
	                	tmp2.Ev1Mv2(h2, o);
	                	double scale = tmp.dot(tmp2);
	                	scale /= tmp.squared();
	                	tmp.TE(scale); //tmp is h1 vector in OH1 direc.
	                	
	                	h2.ME(tmp);
	                	h2.ME(o);
	                	h2.normalize(); // = n normal.
	                	
	                	h2.TE(0.9572*Math.sin(Degree.UNIT.toSim(104.52)));//angle, 104.52, is the SAME!
	                	h2.PEa1Tv1(Math.cos(Degree.UNIT.toSim(104.52)), h1);
	                	h2.PEa1Tv1(1-Math.cos(Degree.UNIT.toSim(104.52)), o);
					}
                	
                	Vector m = a.getPosition();
                	m.E(h1);
                	m.PE(h2);
                	m.PEa1Tv1(-2, o);
                	if(isIce){
                    	m.TE(0.1577/Math.sqrt(m.squared())); //TIP4P/Ice
//                    	m.TE(0.1546/Math.sqrt(m.squared())); //TIP4P2005
                		
                	}else{
                    	m.TE(0.15/Math.sqrt(m.squared())); // TIP4P                		
                	}
                	m.PE(o); 
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
    }


}
