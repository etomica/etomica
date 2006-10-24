package etomica.eigenstuff;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.models.hexane.PairIndexerMolecule;
import etomica.phase.Phase;

public class Reconstituter {
    public static double[][] reconstitute(PairIndexerMolecule pim, Phase phase, double[][][] inData) {
        int dim = phase.space().D();
        //nan let's assume that the leaf atoms are molecules
        double[][] outData = new double[dim*phase.moleculeCount()][dim*phase.moleculeCount()];
        AtomArrayList leafList = phase.getSpeciesMaster().leafList;
        for (int i=0; i<leafList.size(); i++) {
            Atom atom0 = leafList.get(i);
            for (int j=0; j<leafList.size(); j++) {
                Atom atom1 = leafList.get(j);
                int bin = pim.getBin(atom0, atom1);
                for (int k=0; k<dim; k++) {
                    for (int l=0; l<dim; l++) {
                        outData[i*dim+k][j*dim+l] = inData[bin][k][l];
                    }
                }
            }
        }
        return outData;
    }
}
