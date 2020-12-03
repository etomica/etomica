package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesAlkane;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class P4BondTorsionTest {

    Simulation sim;
    PotentialMasterBonding pmBonding;
    int nSpheres;

    @BeforeEach
    void setUp() {
        Space3D space = Space3D.getInstance();
        sim = new Simulation(space);
        nSpheres = 6;
        SpeciesGeneral alkane = SpeciesAlkane.create(nSpheres);
        sim.addSpecies(alkane);
        Box box = sim.makeBox();
        box.setNMolecules(alkane, 1);
        pmBonding = new PotentialMasterBonding(sim, box);
        P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        List<int[]> quads = new ArrayList<>();
        for (int i = 0; i < nSpheres - 3; i++) {
            quads.add(new int[]{i, i + 1, i + 2, i + 3});
        }
        pmBonding.setBondingPotentialQuad(alkane, p4, quads);
        pmBonding.init();
    }

    @Test
    public void testForce() {
        double u = pmBonding.computeAll(true);
        Vector[] f = pmBonding.forces;
        IAtomList atoms = sim.box().getMoleculeList().get(0).getChildList();
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < nSpheres; i++) {
                IAtom a = atoms.get(i);
                Vector3D dr = new Vector3D();
                dr.setRandomSphere(sim.getRandom());
                dr.TE(1e-4);
                Vector fi = new Vector3D();
                fi.E(f[i]);

                a.getPosition().PE(dr);
                double newU = pmBonding.computeAll(true);
                fi.PE(f[i]);
                fi.TE(0.5);
                double duExpected = -fi.dot(dr);
                // we start failing at 1e-10
                assertEquals(duExpected, newU - u, 1e-8);

                u = newU;
            }
        }
    }
}