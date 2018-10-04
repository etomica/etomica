

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.space.Space;

/**
 * Creates a two-body potential group using model parameters.
 *
 * @author navneeth
 */

public class P2PotentialGroupBuilder {

    public static PotentialGroup P2PotentialGroupBuilder(Space space, ModelParams modelParams){

        PotentialGroup potentialGroup = new PotentialGroup(2);
        double sigmaHC = 0.1;

        AtomType[] atomTypes = modelParams.atomTypes;
        double[] sigma = modelParams.sigma;
        double[] epsilon = modelParams.epsilon;
        double[] charge = modelParams.charge;

        for(int i = 0; i < atomTypes.length; i++){
            for(int j = i; j < atomTypes.length; j++){

//                System.out.println("Atom "+i+" and " +j);
//                System.out.println(sigma[i] + " " + sigma[j]);
//                System.out.println(epsilon[i] + " " + epsilon[j]);
//                System.out.println(charge[i] + " " +charge[j]);

                double sigmaij = i == j ? sigma[i]: (sigma[i]+sigma[j])/2;
                double epsilonij= i == j ? epsilon[i]: Math.sqrt(epsilon[i]*epsilon[j]);
                double qiqj = charge[i]*charge[j];
                AtomType[] atomList = new AtomType[] {atomTypes[i],atomTypes[j]};

                P2LennardJones p2LJ;
                if(epsilon[i] != 0 && epsilon[j] != 0) {
                    p2LJ = new P2LennardJones(space, sigmaij, epsilonij);
                    potentialGroup.addPotential(p2LJ, ApiBuilder.makeIntergroupTypeIterator(atomList));
//                    System.out.println("Added p2LJ");
                }

                Potential2SoftSpherical p2ES;
                if(qiqj != 0) {
//                    System.out.print("Added ");
                    if (qiqj < 0 && epsilonij == 0) {
                        p2ES = new P2ElectrostaticWithHardCore(space);
                        ((P2ElectrostaticWithHardCore) p2ES).setCharge1(charge[i]);
                        ((P2ElectrostaticWithHardCore) p2ES).setCharge2(charge[j]);
                        ((P2ElectrostaticWithHardCore) p2ES).setSigma(sigmaHC);
//                        System.out.print("HardCore ");
                    } else {
                        p2ES = new P2Electrostatic(space);
                        ((P2Electrostatic) p2ES).setCharge1(charge[i]);
                        ((P2Electrostatic) p2ES).setCharge2(charge[j]);
                    }
                    potentialGroup.addPotential(p2ES, ApiBuilder.makeIntergroupTypeIterator(atomList));
//                    System.out.println("p2ES");
                }
//                System.out.println();
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
