

package etomica.potential;

import etomica.atom.AtomType;
import etomica.space.Space;
import etomica.species.SpeciesManager;

/**
 * Creates a two-body potential group using model parameters for 1 or 2 species.
 *
 * @author navneeth
 */

public class P2PotentialGroupBuilder {

    public static PotentialMoleculePair P2PotentialGroupBuilder(Space space, SpeciesManager sm, ModelParams MP1, ModelParams MP2){
        boolean debug = false;

        PotentialMoleculePair potentialGroup = new PotentialMoleculePair(space, sm);
        double sigmaHC = 0.1;
        if(MP2 == null) MP2 = MP1;

        for(int i = 0; i < MP1.atomTypes.length; i++){
            int s = MP1 == MP2 ? i : 0;
            for(int j = s; j < MP2.atomTypes.length; j++){

                if(debug) {
                    System.out.println("Atom " + i + " and " + j);
                    System.out.println(MP1.atomTypes[i].getElement() + " " + MP2.atomTypes[j].getElement());
                    System.out.println(MP1.sigma[i] + " " + MP2.sigma[j]);
                    System.out.println(MP1.epsilon[i] + " " + MP2.epsilon[j]);
                    System.out.println(MP1.charge[i] + " " + MP2.charge[j]);
                }
                double sigmaij;
                double epsilonij;
                if( MP1 == MP2 && i == j ){
                    sigmaij = MP1.sigma[i];
                    epsilonij = MP1.epsilon[i];
                }
                else{
                    sigmaij = (MP1.sigma[i]+MP2.sigma[j])/2;
                    epsilonij = Math.sqrt(MP1.epsilon[i]*MP2.epsilon[j]);
                }
                double qiqj = MP1.charge[i]*MP2.charge[j];
                AtomType[] atomList = new AtomType[] {MP1.atomTypes[i],MP2.atomTypes[j]};

                P2LennardJones p2LJ = null;
                if(MP1.epsilon[i] != 0 && MP2.epsilon[j] != 0) {
                    p2LJ = new P2LennardJones(sigmaij, epsilonij);
                    potentialGroup.setAtomPotential(MP1.atomTypes[i], MP2.atomTypes[j], p2LJ);
                    if(debug) {System.out.println("Added p2LJ");}
                }

                IPotential2 p2ES;
                if(qiqj != 0) {
                    if(debug) {System.out.print("Added ");}
                    if (qiqj < 0 && epsilonij == 0) {
                        p2ES = new P2ElectrostaticWithHardCore(MP1.charge[i], MP2.charge[j]);
                        ((P2ElectrostaticWithHardCore) p2ES).setCharge1(MP1.charge[i]);
                        ((P2ElectrostaticWithHardCore) p2ES).setCharge2(MP2.charge[j]);
                        ((P2ElectrostaticWithHardCore) p2ES).setSigma(sigmaHC);
                        if(debug) {System.out.print("HardCore ");}
                    } else {
                        p2ES = new P2Electrostatic(MP1.charge[i], MP2.charge[j]);
                    }
                    IPotential2 p2 = p2ES;
                    if (p2LJ != null) p2 = new P2SoftSphericalSum(p2LJ, p2ES);
                    potentialGroup.setAtomPotential(MP1.atomTypes[i], MP2.atomTypes[j], p2);
                    if(debug) {System.out.println("p2ES");}
                }
                if(debug) {System.out.println();}
            }
        }
        return potentialGroup;
    }

    public static class ModelParams {

        protected AtomType[] atomTypes;
        protected double[] sigma;
        protected double[] epsilon;
        protected double[] charge;

        public ModelParams(AtomType[] atomTypes, double[] sigma, double[] epsilon, double[] charge){
            this.atomTypes = atomTypes;
            this.sigma = sigma;
            this.epsilon = epsilon;
            this.charge = charge;
        }
    }

}
