package etomica.virial.simulations.mbnrg;

import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;

import java.util.ArrayList;

public class energy1b implements IPotentialMolecular {
    /**
     * @brief Gets the one body energy for a given set of monomers of the same
     * monomer type.
     *
     * Given the monomer id and the number of monomers, will return the total sum
     * of the 1b energy of those monomers.
     * @param[in] mon Monomer id
     * @param[in] nm number of monomers of monomer type "mon"
     * @param[in] xyz coordinates of the monomer
     * unittest-energy1b.cpp has sample xyz1
     * @param[in,out] bad_idxs Vector with the indexes o extremely distorted monomers
     * has an energy larger than the value set in definitions.h (EMAX1B)
     * @return Sum of the one-body energies of all the monomers passed as arguments
     */
    static double get_1b_energy(String mon1, Integer nm, Double[] xyz1, ArrayList<Integer> bad_idxs) {

        ArrayList<Double> energies;
        // Look for the proper call to energy depending on the monomer id
        if (mon1.equals("co2")) {
//            x1b_A1B2_deg4::x1b_A1B2_v1x pot(mon1);

            A1B2_1b pot = new A1B2_1b("co2");
            energies = pot.eval(xyz1, nm);
        }
        else {
            return 0.0;
        }

        // Add total energy, and check for too high energies
        double e = 0.0;
        for (int i = 0; i < nm; i++) {
            e += energies.get(i);
            if (energies.get(i) > Definitions.EMAX1B) bad_idxs.add(i);
        }

        // Return energy
        return e;
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return 0;
    }
//    static double get_1b_energy(String mon1, Integer nm, ArrayList<Double> xyz1, ArrayList<Double> grad1,
//                         ArrayList<Integer> bad_idxs) {
//
//        ArrayList<Double> energies;
//        // Look for the proper call to energy depending on the monomer id
//        if (mon1 == "co2") {
//            A1B2 pot = new A1B2();
//            energies = pot.eval(xyz1, nm);
//        } else {
//            return 0.0;
//        }
//
//        // Add total energy, and check for too high energies
//        double e = 0.0;
//        for (int i = 0; i < nm; i++) {
//            e += energies.get(i);
//            if (energies.get(i) > Definitions.EMAX1B) bad_idxs.add(i);
//        }
//
//        // Return energy
//        return e;
//    }

}
