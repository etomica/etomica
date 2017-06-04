/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.space.Boundary;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Degree;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Determines c/a ratio of primitive vector lengths corresponding to the
 * minimum energy HCP structure.
 * 
 * @author Andrew Schultz
 */
public class MinimizeHCP extends Simulation {

    private static final long serialVersionUID = 1L;
    public final Box box;
    public final PotentialMaster potentialMaster;
    public final double density;
    public final int numMolecules;
    public MinimizeHCP(Space _space, int numAtoms, double density, int exponent, double rc, double initC) {
        super(_space);


        this.density = density;
        this.numMolecules = numAtoms;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        setC(initC, density);

        if(exponent ==6){
        	System.out.println("set rc");
        	rc = box.getBoundary().getBoxSize().getX(0)*0.495;
        }

        BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, space);
        BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class);
        potentialMaster = new PotentialMasterList(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, space), space);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        potential = new P2SoftSphericalTruncated(space, potential, rc);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(false);


        ((PotentialMasterList)potentialMaster).setRange(rc);
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        ((PotentialMasterList)potentialMaster).getNeighborManager(box).reset();

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*box.getBoundary().getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see MinimizeHCP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        int exponentN = params.exponentN;
        final int numMolecules = params.numMolecules;
        double initC = params.initC;
        double rc = params.rc;
        double mincf = params.mincf;
        double maxcf = params.maxcf;

        System.out.println("Running HCP soft sphere overlap simulation");
        System.out.println(numMolecules + " atoms at density " + density);
        System.out.println("exponent N: " + exponentN);

        MinimizeHCP sim = new MinimizeHCP(Space.getInstance(3), numMolecules, density, exponentN, rc, initC);

        double cfvalue = sim.findcf(initC, maxcf, mincf);
        System.out.println("final cf: " + cfvalue);

    }

    public void setC(double newC, double density) {
        double a = Math.pow(Math.sqrt(2)/density, 1.0/3.0);
        double c = newC*a;
        a *=  Math.sqrt(Math.sqrt(8.0/3.0)/newC);
        //System.out.println("cf "+newC);
        int nC = (int)Math.ceil(Math.pow(box.getLeafList().getAtomCount()/2, 1.0/3.0));
        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{nC*a, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-nC*a*Math.cos(Degree.UNIT.toSim(60)), nC*a*Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, nC*c});

        Primitive primitive = new PrimitiveHexagonal(space, nC*a, nC*c);

        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        Basis basisHCP = new BasisHcp();

        box.setBoundary(boundary);

        ConfigurationLattice config = new ConfigurationLattice(new BravaisLatticeCrystal(primitive, basisHCP), space);
        config.initializeCoordinates(box);
    }

    public double findcf(double initC, double maxcf, double mincf){

        int bootstrap = 0;
    	double[] u = new double[3];
    	double[] allcf = new double[3];

        double cf = mincf;

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
    	meterPE.setBox(box);


        while (true) {

            setC(cf, density);

            double latticeEnergy = meterPE.getDataAsScalar();
//            System.out.println("lattice energy: " + latticeEnergy/numMolecules);
            if (bootstrap < 3) {
                allcf[bootstrap] = cf;
                u[bootstrap] = latticeEnergy;
                bootstrap++;
                cf += 0.5*(maxcf-mincf);
            }
            else {
                if (cf > allcf[2]) {
                    allcf[0] = allcf[1];
                    allcf[1] = allcf[2];
                    allcf[2] = cf;
                    u[0] = u[1];
                    u[1] = u[2];
                    u[2] = latticeEnergy;
                }
                else if (cf < allcf[0]) {
                    allcf[2] = allcf[1];
                    allcf[1] = allcf[0];
                    allcf[0] = cf;
                    u[2] = u[1];
                    u[1] = u[0];
                    u[0] = latticeEnergy;
                }
                else if (u[2] > u[0]) {
                    maxcf = allcf[2];
                    if (cf > allcf[1]) {
                        u[2] = latticeEnergy;
                        allcf[2] = cf;
                    }
                    else {
                        u[2] = u[1];
                        allcf[2] = allcf[1];
                        u[1] = latticeEnergy;
                        allcf[1] = cf;
                    }
                }
                else {
                    mincf = allcf[0];
                    if (cf < allcf[1]) {
                        u[0] = latticeEnergy;
                        allcf[0] = cf;
                    }
                    else {
                        u[0] = u[1];
                        allcf[0] = allcf[1];
                        u[1] = latticeEnergy;
                        allcf[1] = cf;
                    }
                }

//                System.out.println("all cf:");
//                System.out.println(allcf[0]+" "+u[0]/numMolecules);
//                System.out.println(allcf[1]+" "+u[1]/numMolecules);
//                System.out.println(allcf[2]+" "+u[2]/numMolecules);
                if (u[1] > u[0] && u[1] > u[2]) {
                    // we found a maximum, due to numerical precision failure
                    // just bail and pretend that the middle point is the global minimum
                  //  System.out.println("oops, a maximum");
//                  System.out.println("all cf:");
//                  System.out.println(allcf[0]+" "+u[0]/numMolecules);
//                  System.out.println(allcf[1]+" "+u[1]/numMolecules);
//                  System.out.println(allcf[2]+" "+u[2]/numMolecules);
//                    System.out.println("c/a "+allcf[1]);
                    return allcf[1];
                }
            }

            if (bootstrap == 3) {
                // now estimate minimum in U from the three points.
                double dc01 = allcf[1]-allcf[0];
                double dc12 = allcf[2]-allcf[1];
                double du01 = u[1]-u[0];
                double du12 = u[2]-u[1];
                double dudc01 = du01/dc01;
                double dudc12 = du12/dc12;
                double m = (dudc12-dudc01)/(0.5*(dc01+dc12));
                cf = 0.9*(0.5*(allcf[1]+allcf[2]) - dudc12/m) + 0.1*(0.5*(allcf[0]+allcf[2]));
                if (cf == allcf[1] || cf == allcf[2]) {
                    cf = 0.5*(allcf[1] + allcf[2]);
                }
                if (cf == allcf[0] || cf == allcf[1]) {
                    cf = 0.5*(allcf[1] + allcf[0]);
                }
                if (cf < mincf) {
                    cf = 0.5*(mincf + allcf[0]);
                }
                if (cf > maxcf) {
                    cf = 0.5*(maxcf + allcf[2]);
                }

                if (cf == allcf[0] || cf == allcf[1] || cf == allcf[2]) {
                    // we converged cf to numerical precision.
                    //System.out.println("c/a "+cf);
                    return cf;
                }
            }
        }
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 3456;
        public double density = 2.20617;
        public int exponentN = 6;
        public double rc = 2.2;
        public double initC = 1.6339554555247535;//Math.sqrt(8.0/3.0);
        public double mincf = 1.2;
        public double maxcf = 2.0;
    }
}
