package etomica.modules.reactionequilibrium;

import etomica.Atom;
import etomica.MeterAbstract;
import etomica.Phase;
import etomica.Simulation;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.units.Dimension;

public final class MeterDimerFraction extends MeterAbstract {
    String[] labels = new String[5];
    double[] currentValues = new double[5];
    int[] count = new int[5];
    public final int idx;
    private AtomIteratorListSimple iterator = new AtomIteratorListSimple();
    public MeterDimerFraction(Simulation sim, int idx) {
        super(5);
        labels[0] = "R";
        labels[1] = "B";
        labels[2] = "R-R";
        labels[3] = "R-B";
        labels[4] = "B-B";
        this.idx = idx;
    }
    
    public double[] getData(Phase phase) {
        for(int i=0; i<count.length; i++) {count[i] = 0;}
        iterator.setList(phase.speciesMaster.atomList);
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
        for(int i=0; i<count.length; i++) {
        	currentValues[i] = count[i]/nMole;
        }
        currentValues[2] *= 0.5;
        currentValues[3] *= 0.5;
        currentValues[4] *= 0.5;
        return currentValues;
    }
    
    public Dimension getXDimension() {return Dimension.NULL;}
    public Dimension getDimension() {return Dimension.NULL;}
}