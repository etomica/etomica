/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.reactionequilibrium;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.ISpecies;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Fraction;
import etomica.units.Quantity;
import etomica.units.Volume;

public final class MeterDimerFraction implements IEtomicaDataSource {
    public MeterDimerFraction(AtomLeafAgentManager aam) {
        data = new DataTable(1,5);
        DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Dimer Fraction", Fraction.DIMENSION, new int[]{5});
        dataInfo = new DataInfoTable("Dimer Fraction", new DataInfoDoubleArray[]{columnInfo}, 5, new String[]{"R", "B", "R-R", "R-B", "B-B"});
        agentManager = aam;
        tag = new DataTag();
        dataInfo.addTag(tag);
        dataDensity = new DataDouble();
        dataGroup = new DataGroup(new IData[] {data, dataDensity});
        DataDouble.DataInfoDouble dataInfoDouble = new DataInfoDouble("rho", new CompoundDimension(new Dimension[] {Quantity.DIMENSION, Volume.DIMENSION},new double[] {1,-1}));
        dataInfoGroup = new DataInfoGroup("x and rho", Dimension.MIXED, new IEtomicaDataInfo[] {dataInfo, dataInfoDouble});
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfoGroup;
    }

    public DataTag getTag() {
        return tag;
    }
    
    public void setSpeciesA(ISpecies newSpeciesA) {
        speciesA = newSpeciesA;
    }
    
    public IData getData() {
        for(int i=0; i<count.length; i++) {count[i] = 0;}
        IAtomList leafAtoms = box.getLeafList();
        int nLeaf = leafAtoms.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            IAtom a = leafAtoms.getAtom(i);
        	IAtom partner = (IAtom)agentManager.getAgent(a);
  //      	if(partner != null) System.out.println(a.node.index()+" "+partner.node.index());
            if(a.getType().getSpecies() == speciesA) {
               if(partner == null) {
                 count[0]++;  //A radical
               }
               else if(partner.getType().getSpecies() == speciesA) {
                 count[2]++;  //A-A
               }
               else {
                 count[3]++;  //A-B
               }
            }
            else { //a is species 2
               if(partner == null) {
                 count[1]++;  //B radical
               }
               else if(partner.getType().getSpecies() == speciesA) {
                 count[3]++;  //A-B
               }
               else {
                 count[4]++;  //B-B
               }
            }//end of if
        
        }//end of for loop
        
        double nMole = count[0] + count[1] + 0.5*(count[2]+count[3]+count[4]);
        double[] x = ((DataDoubleArray)data.getData(0)).getData();
        for(int i=0; i<count.length; i++) {
        	x[i] = count[i]/nMole;
        }
        x[2] *= 0.5;
        x[3] *= 0.5;
        x[4] *= 0.5;
        
        dataDensity.x = nMole/box.getBoundary().volume();
        return dataGroup;
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
        iterator.setBox(box);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private String name;
    private Box box;
    private ISpecies speciesA;
    private final DataTable data;
    private final DataInfoTable dataInfo;
    private int[] count = new int[5];
    private AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    protected final AtomLeafAgentManager agentManager;
    private final DataTag tag;
    private final DataGroup dataGroup;
    private final DataDouble dataDensity;
    private final DataGroup.DataInfoGroup dataInfoGroup;
}
