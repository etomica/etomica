/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.integrator.mcmove.MCMoveVolumeFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.IFunction;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3DFasterer extends Simulation {

    public final PotentialMasterListFasterer potentialMasterList;
    public final PotentialMasterFasterer potentialMasterLongCut;
    public IntegratorVelocityVerletFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public Potential2SoftSpherical potential;
    public IntegratorMCFasterer integratorMC;
    public MCMoveVolumeFasterer mcMoveVolume;


    public LjMd3DFasterer(int numAtoms, double temperature, double density, double pressure, double tStep, double rcShort, double rcLong, int hybridInterval, IFunction vBias, boolean ss) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        double L = Math.pow(numAtoms / density, 1.0 / 3.0);
        double nbrRange = rcShort * 1.6;
        if (nbrRange > 0.5 * L) {
            if (rcShort > 0.4 * L) {
                throw new RuntimeException("rcShort is too large");
            }
            nbrRange = 0.495 * L;
        }
        potentialMasterList = new PotentialMasterListFasterer(getSpeciesManager(), box, 2, nbrRange, BondingInfo.noBonding());
        potentialMasterList.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerletFasterer(potentialMasterList, random, tStep, temperature, box);
        integrator.setIsothermal(true);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        potential = ss ? new P2SoftSphere(space, 1, 4, 12) : new P2LennardJones(space);
        AtomType leafType = species.getLeafType();
        TruncationFactory f;
        P2SoftSphericalTruncatedForceShifted potentialTruncatedForceShifted = new P2SoftSphericalTruncatedForceShifted(space, potential, rcShort);

        potentialMasterList.setPairPotential(leafType, leafType, potentialTruncatedForceShifted);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.HYBRID_MC);
        integrator.setThermostatInterval(hybridInterval);

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);

        double rc = rcLong > 0 ? rcLong : 0.495 * box.getBoundary().getBoxSize().getX(0);
        while (rc >= 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            rc -= 0.5;
            System.out.println("long rc => " + rc);
            //throw new RuntimeException("rc must be less than half the box");
        }
        if (rcLong <= 0) {
            System.out.println("long rc: " + rc);
        }

        // untruncated.  users of this will truncate on their own
        potentialMasterLongCut = new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());

        potentialMasterLongCut.setPairPotential(leafType, leafType, potential);

        if (!Double.isNaN(pressure)) {
            // XXX neighbor list and MC?  don't we miss stuff when we compress?
            integratorMC = new IntegratorMCFasterer(potentialMasterList, random, temperature, box);
            mcMoveVolume = new MCMoveVolumeFasterer(integratorMC, random, pressure);
            mcMoveVolume.setVolumeBias(vBias);
            integratorMC.getMoveManager().addMCMove(mcMoveVolume);
            integrator.setIntegratorMC(integratorMC, 1);
        }
    }
}
