/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.thermalconductivity;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.IPotential1;
import etomica.potential.IPotential2;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Heat conduction through a 2D fluid of hard disks (or an ideal gas) held between a hot
 * wall (x = 0) and a cold wall (x = Lx). The walls are thermal (diffuse) walls, the box
 * is periodic in y, and there is no thermostat in the interior. Reduced units
 * (sigma = m = kB = 1).
 */
public class ThermalConductivity extends Simulation {

    public static final double SIGMA = 1.0;

    public final SpeciesGeneral species;
    public final Box box;
    public final IntegratorHard integrator;
    public final P2HardGeneric p2hard;
    public final P1HardThermalWall wall;
    public final ConductionMeter meter;

    protected final double lx, ly;
    protected double phi;
    protected boolean idealGas = false;

    public ThermalConductivity() {
        this(40, 20, 0.2, 2.0, 1.0);
    }

    public ThermalConductivity(double lx, double ly, double phi, double tHot, double tCold) {
        super(Space2D.getInstance());
        this.lx = lx;
        this.ly = ly;
        this.phi = phi;

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) species.getLeafType().getElement()).setMass(1.0);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularSlit(0, space));
        box.getBoundary().setBoxSize(Vector.of(lx, ly));

        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        int nTypes = getSpeciesManager().getAtomTypeCount();
        IPotential2[][] pairPotentials = new IPotential2[nTypes][nTypes];
        IPotential1[] fieldPotentials = new IPotential1[nTypes];

        p2hard = P2HardSphere.makePotential(SIGMA);
        pairPotentials[species.getLeafType().getIndex()][species.getLeafType().getIndex()] = p2hard;
        wall = new P1HardThermalWall(box, random, tHot, tCold, 0.5 * SIGMA);
        fieldPotentials[species.getLeafType().getIndex()] = wall;

        integrator = new IntegratorHard(pairPotentials, fieldPotentials, neighborManager, random, 0.05,
                0.5 * (tHot + tCold), box, getSpeciesManager(), null);
        // no thermostat: the walls are the only heat source and sink
        integrator.setIsothermal(false);
        integrator.setMaxCollisionDiameter(species.getLeafType(), SIGMA);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        getController().addActivity(new ActivityIntegrate(integrator));

        meter = new ConductionMeter(this);
        integrator.getEventManager().addListener(meter);

        box.setNMolecules(species, numberForPhi(phi));
        initConfiguration();
    }

    public double getLx() {
        return lx;
    }

    public double getLy() {
        return ly;
    }

    public double getPhi() {
        return phi;
    }

    public boolean isIdealGas() {
        return idealGas;
    }

    public double getMeanTemperature() {
        return 0.5 * (wall.getTHot() + wall.getTCold());
    }

    public int numberForPhi(double phi) {
        return Math.max(1, (int) Math.round(4 * phi * lx * ly / (Math.PI * SIGMA * SIGMA)));
    }

    /**
     * Sets the area fraction (changing N) and reinitializes the configuration.
     * Must be called from the controller thread.
     */
    public void setPhi(double newPhi) {
        phi = newPhi;
        box.setNMolecules(species, numberForPhi(newPhi));
        initConfiguration();
        integrator.reset();
        meter.reset();
    }

    /**
     * Selects ideal gas (disks pass through each other) or hard disks and reinitializes
     * the configuration.  Must be called from the controller thread.
     */
    public void setIdealGas(boolean ideal) {
        idealGas = ideal;
        AtomType type = species.getLeafType();
        integrator.setPairPotential(type, type, ideal ? null : p2hard);
        initConfiguration();
        integrator.reset();
        meter.reset();
    }

    /**
     * Sets the wall temperatures.  Must be called from the controller thread.
     */
    public void setWallTemperatures(double tHot, double tCold) {
        wall.setTemperatures(tHot, tCold);
        meter.reset();
    }

    /**
     * Configuration used by the Reinitialize button.
     */
    public Configuration makeConfiguration() {
        return b -> {
            initConfiguration();
            meter.reset();
        };
    }

    /**
     * Places the disks on a random subset of the sites of a rectangular grid that
     * keeps all disks off the walls and (for hard disks) apart from each other, and
     * gives them Maxwell-Boltzmann velocities at the mean wall temperature.
     */
    public void initConfiguration() {
        IAtomList atoms = box.getLeafList();
        int n = atoms.size();
        if (n == 0) return;
        double minSpacing = 1.02 * SIGMA;
        int nxMax = Math.max(1, (int) Math.floor((lx - SIGMA) / minSpacing));
        int nyMax = Math.max(1, (int) Math.floor(ly / minSpacing));
        int ny = (int) Math.round(Math.sqrt(n * ly / (lx - SIGMA)));
        ny = Math.min(Math.max(ny, 1), nyMax);
        int nx = (int) Math.ceil((double) n / ny);
        if (nx > nxMax) {
            nx = nxMax;
            ny = (int) Math.ceil((double) n / nx);
        }
        double dx = (lx - SIGMA) / nx;
        double dy = ly / ny;
        int nSites = nx * ny;
        int[] sites = new int[nSites];
        for (int i = 0; i < nSites; i++) sites[i] = i;
        // shuffle so the (possibly partially filled) grid is uniformly occupied
        for (int i = nSites - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int tmp = sites[i];
            sites[i] = sites[j];
            sites[j] = tmp;
        }
        double sd = Math.sqrt(getMeanTemperature() / species.getLeafType().getMass());
        double[] vSum = new double[2];
        for (int i = 0; i < n; i++) {
            IAtomKinetic atom = (IAtomKinetic) atoms.get(i);
            // if there are more atoms than sites (ideal gas only), recycle sites
            int s = sites[i % nSites];
            int ix = s / ny;
            int iy = s % ny;
            atom.getPosition().setX(0, -0.5 * lx + 0.5 * SIGMA + (ix + 0.5) * dx);
            atom.getPosition().setX(1, -0.5 * ly + (iy + 0.5) * dy);
            for (int k = 0; k < 2; k++) {
                double v = sd * random.nextGaussian();
                atom.getVelocity().setX(k, v);
                vSum[k] += v;
            }
        }
        // remove net drift
        for (int i = 0; i < n; i++) {
            IAtomKinetic atom = (IAtomKinetic) atoms.get(i);
            for (int k = 0; k < 2; k++) {
                atom.getVelocity().setX(k, atom.getVelocity().getX(k) - vSum[k] / n);
            }
        }
    }

    public static void main(String[] args) {
        ThermalConductivity sim = new ThermalConductivity();
        double tEquil = 400, tRun = 4000;
        if (args.length > 0) sim.setPhi(Double.parseDouble(args[0]));
        sim.integrator.reset();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, (long) (tEquil / 0.05)));
        sim.meter.reset();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, (long) (tRun / 0.05)));
        sim.meter.printSummary();
    }
}
