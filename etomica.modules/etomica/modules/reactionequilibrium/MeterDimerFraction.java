package etomica.modules.reactionequilibrium;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.phase.Phase;
import etomica.units.Dimension;
import etomica.utility.NameMaker;

public final class MeterDimerFraction implements DataSource, Meter {
    public MeterDimerFraction(int idx) {
        data = new DataDoubleArray("Dimer Fraction",Dimension.FRACTION,5);
        setName(NameMaker.makeName(this.getClass()));
        labels[0] = "R";
        labels[1] = "B";
        labels[2] = "R-R";
        labels[3] = "R-B";
        labels[4] = "B-B";
        this.idx = idx;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

    public Data getData() {
        for(int i=0; i<count.length; i++) {count[i] = 0;}
        iterator.reset();
        while(iterator.hasNext()) {
        	Atom a = iterator.nextAtom();
        	Atom partner = (Atom)a.allatomAgents[idx];
  //      	if(partner != null) System.out.println(a.node.index()+" "+partner.node.index());
            if(a.type.getSpeciesIndex()== 0) {
               if(partner == null) {
                 count[0]++;  //A radical
               }
               else if(partner.type.getSpeciesIndex()== 0) {
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
               else if(partner.type.getSpeciesIndex()== 0) {
                 count[3]++;  //A-B
               }
               else {
                 count[4]++;  //B-B
               }
            }//end of if
        
        }//end of for loop
        
        double nMole = count[0] + count[1] + 0.5*(count[2]+count[3]+count[4]);
        double[] x = data.getData();
        for(int i=0; i<count.length; i++) {
        	x[i] = count[i]/nMole;
        }
        x[2] *= 0.5;
        x[3] *= 0.5;
        x[4] *= 0.5;
        return data;
    }
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        iterator.setPhase(phase);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private String name;
    private Phase phase;
    private final DataDoubleArray data;
    private String[] labels = new String[5];
    private int[] count = new int[5];
    public final int idx;
    private AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();

}