package etomica.modules.reactionequilibrium;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.phase.Phase;
import etomica.units.Fraction;
import etomica.util.NameMaker;

public final class MeterDimerFraction implements Meter {
    public MeterDimerFraction(ReactionEquilibrium sim) {
        data = new DataTable(1,5);
        DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Dimer Fraction", Fraction.DIMENSION, new int[]{5});
        dataInfo = new DataInfoTable("Dimer Fraction", new DataInfoDoubleArray[]{columnInfo}, 5, new String[]{"R", "B", "R-R", "R-B", "B-B"});
        setName(NameMaker.makeName(this.getClass()));
        agentSource = sim;
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        agents = agentSource.getAgents(phase);
        for(int i=0; i<count.length; i++) {count[i] = 0;}
        iterator.reset();
        while(iterator.hasNext()) {
        	Atom a = iterator.nextAtom();
        	Atom partner = agents[a.getGlobalIndex()];
  //      	if(partner != null) System.out.println(a.node.index()+" "+partner.node.index());
            if(a.getType().getSpeciesIndex()== 1) {
               if(partner == null) {
                 count[0]++;  //A radical
               }
               else if(partner.getType().getSpeciesIndex()== 1) {
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
               else if(partner.getType().getSpeciesIndex()== 1) {
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
    private final DataTable data;
    private final DataInfoTable dataInfo;
    private int[] count = new int[5];
    private AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    protected final ReactionEquilibrium agentSource;
    protected Atom[] agents;
    private final DataTag tag;
}