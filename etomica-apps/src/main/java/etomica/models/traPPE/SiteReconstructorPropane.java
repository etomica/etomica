/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.traPPE;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationLinear;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMoleculePairImplicit;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.util.Arrays;

/**
 * Site reconstructor for TraPPE-EH propane.
 *
 * Child-list order for each molecule:
 *   0: C1   (left methyl carbon)
 *   1: C2   (methylene carbon)
 *   2: C3   (right methyl carbon)
 *   3: HL   (one hydrogen on C1)
 *   4: HR   (one hydrogen on C3)
 *
 * From these, this class reconstructs all implied LJ interaction sites:
 *   explicit carbon sites: C1, C2, C3
 *   left methyl bond-center sites: 3
 *   methylene bond-center sites:   2
 *   right methyl bond-center sites:3
 *
 * Total: 11 LJ sites per molecule.
 */
public class SiteReconstructorPropane implements PotentialMoleculePairImplicit.SiteReconstructor {

    // Geometry constants
    protected static final double DCH = 1.100;  // Angstrom
    protected static final double HALF_DCH = 0.5 * DCH;
    protected static final double COS120 = -0.5;
    protected static final double SIN120 = Math.sqrt(3)/2;
    public final ISpecies species;
    public SiteReconstructorPropane(ISpecies species){
        this.species = species;
    }
    public PotentialMoleculePairImplicit.SiteSet makeSites() {
        AtomType[] types = new AtomType[11];
        // first 3 sites are C, then 8 CH (H are not included)
        Arrays.fill(types, 0, 3, species.getAtomType(0));
        Arrays.fill(types, 3, 11, species.getAtomType(2));
        return new PotentialMoleculePairImplicit.SiteSet(types);
    }


    public void reconstructSites(IMolecule molecule, PotentialMoleculePairImplicit.SiteSet sites) {
        IAtomList atoms = molecule.getChildList();

        IAtom aC1 = atoms.get(0);
        IAtom aC2 = atoms.get(1);
        IAtom aC3 = atoms.get(2);
        IAtom aHL = atoms.get(3);
        IAtom aHR = atoms.get(4);

        Vector rC1 = aC1.getPosition();
        Vector rC2 = aC2.getPosition();
        Vector rC3 = aC3.getPosition();
        Vector rHL = aHL.getPosition();
        Vector rHR = aHR.getPosition();

        // 0-2: carbon LJ sites
        sites.pos[0].E(rC1);
        sites.pos[1].E(rC2);
        sites.pos[2].E(rC3);

        // 3-5: left methyl bond-center sites
        reconstructMethylBondCenters(
                rC1, rC2, rHL,
                sites.pos[3], sites.pos[4], sites.pos[5]
        );

        // 6-7: central methylene bond-center sites
        reconstructMethyleneBondCenters(
                rC1, rC2, rC3,
                sites.pos[6], sites.pos[7]
        );

        // 8-10: right methyl bond-center sites
        reconstructMethylBondCenters(
                rC3, rC2, rHR,
                sites.pos[8], sites.pos[9], sites.pos[10]
        );
        for (Vector v:sites.pos){
            if (v.squared() > 10000) {
                System.out.println("faraway " + molecule);
                for (Vector v2:sites.pos){
                    System.out.println("v " + v2);
                }
                for (IAtom a:atoms) {
                    System.out.println("a " + a.getLeafIndex() + " " + a.getPosition());}
                PotentialMoleculePairImplicit.debug = true;
                return;
            }
        }
    }

    /**
     * Reconstruct the three bond-center LJ sites on a terminal methyl carbon.
     *
     * terminalCarbon = methyl carbon
     * neighborCarbon = adjacent backbone carbon
     * knownH         = one explicit hydrogen on that methyl
     *
     * Output bc1, bc2, bc3 are the three C-H bond-center LJ site positions.
     */
    protected static void reconstructMethylBondCenters(
            Vector terminalCarbon,
            Vector neighborCarbon,
            Vector knownH,
            Vector bc1,
            Vector bc2,
            Vector bc3) {

        Vector3D axis = new Vector3D();
        axis.Ev1Mv2(neighborCarbon, terminalCarbon);
        axis.normalize();

        Vector3D u1 = new Vector3D();
        u1.Ev1Mv2(knownH, terminalCarbon);
        u1.normalize();

        Vector3D u2 = rotate120(u1, axis, +1.0);
        Vector3D u3 = rotate120(u1, axis, -1.0);

        // Bond-center sites are at carbon + 0.5*dCH*u
        bc1.E(terminalCarbon);
        bc1.PEa1Tv1(HALF_DCH, u1);

        bc2.E(terminalCarbon);
        bc2.PEa1Tv1(HALF_DCH, u2);

        bc3.E(terminalCarbon);
        bc3.PEa1Tv1(HALF_DCH, u3);
    }

    /**
     * Reconstruct the two bond-center LJ sites on the central methylene carbon.
     *
     * Uses the symmetric CH2 construction:
     *   - the two H directions are symmetric about the direction opposite the
     *     bisector of the two C-C bonds
     *   - they are displaced equally above/below the C1-C2-C3 plane
     *
     * This gives the two central C-H bond-center LJ sites.
     */
    protected static final double THETA_HCH_CH2 = Math.toRadians(107.8);
    protected static final double COS_HALF_HCH = Math.cos(0.5 * THETA_HCH_CH2);
    protected static final double SIN_HALF_HCH = Math.sin(0.5 * THETA_HCH_CH2);

    protected static void reconstructMethyleneBondCenters(
            Vector rC1,
            Vector rC2,
            Vector rC3,
            Vector bcA,
            Vector bcB) {

        // b1 = unit vector from C2 toward C1
        Vector3D b1 = new Vector3D();
        b1.Ev1Mv2(rC1, rC2);
        b1.normalize();

        // b3 = unit vector from C2 toward C3
        Vector3D b3 = new Vector3D();
        b3.Ev1Mv2(rC3, rC2);
        b3.normalize();

        // s = bisector in the C-C-C plane, pointing away from the two carbons
        Vector3D s = new Vector3D();
        s.Ev1Pv2(b1, b3);
        if (s.squared() < 1.0e-14) {
            throw new RuntimeException("Degenerate methylene geometry: bisector nearly zero");
        }
        s.normalize();
        s.TE(-1.0);  // point away from the two carbons

        // n = normal to the C1-C2-C3 plane
        Vector3D n = new Vector3D();
        n.E(b1);
        n.XE(b3);
        if (n.squared() < 1.0e-14) {
            // fallback if nearly collinear
            n.setPerpendicularTo(s);
        }
        n.normalize();

        // uA = cos(alpha) s + sin(alpha) n
        Vector3D uA = new Vector3D();
        uA.Ea1Tv1(COS_HALF_HCH, s);
        uA.PEa1Tv1(SIN_HALF_HCH, n);

        // uB = cos(alpha) s - sin(alpha) n
        Vector3D uB = new Vector3D();
        uB.Ea1Tv1(COS_HALF_HCH, s);
        uB.PEa1Tv1(-SIN_HALF_HCH, n);

        // Bond-center sites = C2 + 0.5*dCH*u
        bcA.E(rC2);
        bcA.PEa1Tv1(HALF_DCH, uA);

        bcB.E(rC2);
        bcB.PEa1Tv1(HALF_DCH, uB);
    }

    /**
     * Rotate u by +/-120 degrees about the unit axis.
     * sign = +1 => +120 degrees
     * sign = -1 => -120 degrees
     *
     * Rodrigues formula:
     *   u' = u cos t + (axis x u) sin t + axis (axis·u)(1-cos t)
     */
    protected static Vector3D rotate120(Vector u, Vector axis, double sign) {
        Vector3D term1 = new Vector3D();
        term1.Ea1Tv1(COS120, u);

        Vector3D cross = new Vector3D();
        cross.E(axis);
        cross.XE(u);

        Vector3D term2 = new Vector3D();
        term2.Ea1Tv1(sign * SIN120, cross);

        Vector3D term3 = new Vector3D();
        term3.E(axis);
        term3.TE(axis.dot(u) * (1.0 - COS120));

        Vector3D out = new Vector3D();
        out.E(term1);
        out.PE(term2);
        out.PE(term3);

        return out;
    }

    public static void main(String[] args) {
        Space3D s = Space3D.getInstance();
        SpeciesBuilder sb = new SpeciesBuilder(s);
        AtomType C = new AtomType(Carbon.INSTANCE);
        AtomType H = new AtomType(Hydrogen.INSTANCE);
        AtomType CH = new AtomType(new ElementSimple("CH", 0));
        sb.addCount(C, 3);
        sb.addCount(H, 2);
        sb.addCount(CH, 8);
        sb.withConformation(new ConformationLinear(s));
        Box box = new Box(s);
        Simulation sim = new Simulation(s);
        ISpecies species = sb.build();
        sim.addSpecies(species);
        sim.addBox(box);
        box.setNMolecules(species, 1);

        int[][] bonds = new int[][]{{0,1},{1,2},
                {0,3},{0,4},{0,5},
                {1,6},{1,7},
                {2,8},{2,9},{2,10}};
        int[][] angles = new int[][]{{0,1,2},  // C0C1C2, free
                {3,0,1},{4,0,1},{5,0,1},       // HC0C1, should be 110.7
                {0,1,6},{0,1,7},{2,1,6},{2,1,7}, // C0C1H and C2C1H, value depends on CCC, should all be equal
                {1,2,8},{1,2,9},{1,2,10},      // C1C2H, should be 110.7
                {3,0,4},{3,0,5},{4,0,5},       // HC0H, should be 107.8
                {8,2,9},{8,2,10},{9,2,10},     // HC2H, should be 107.8
                {6,1,7}                        // HC1H, should be 107.8
        };
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        int[][] bonds0 = new int[][]{{0,1},{1,2},{0,3},{2,4}};
        double[] bl0 = new double[]{1.535,1.535,1.1,1.1};
        int[][] angles0 = new int[][]{{1,0,3},{1,2,4}};
        double[] phi0 = new double[]{Math.PI*110.70/180, Math.PI*110.7/180};
        PotentialMoleculePairImplicit.SiteReconstructor reconstructor = new SiteReconstructorPropane(species);

        for (int step=0; step<10; step++) {
            for (int i=0; i<bonds0.length; i++) {
                int[] b0 = bonds0[i];
                Vector r = atoms.get(b0[1]).getPosition();
                r.setRandomSphere(sim.getRandom());
                for (int j=0; j<angles0.length; j++) {
                    // look for a bond angle that includes this atom
                    if (angles0[j][2] == b0[1]) {
                        // we have a bond angle, rotate r to match the prescribed bond angle
                        Vector x = s.makeVector();
                        x.Ev1Mv2(atoms.get(angles0[j][0]).getPosition(), atoms.get(angles0[j][1]).getPosition());
                        x.normalize();
                        Vector y = s.makeVector();
                        y.E(r);
                        double xydot = x.dot(y);
                        y.PEa1Tv1(-xydot, x);
                        y.normalize();
                        r.Ea1Tv1(Math.cos(phi0[j]), x);
                        r.PEa1Tv1(Math.sin(phi0[j]), y);
                        break;
                    }
                }
                // now set our bond length and shift
                r.TE(bl0[i]);
                r.PE(atoms.get(b0[0]).getPosition());
            }
            PotentialMoleculePairImplicit.SiteSet sites = reconstructor.makeSites();

            reconstructor.reconstructSites(box.getMoleculeList().get(0), sites);
            System.out.println();
            for (int j = 0; j < bonds.length; j++) {
                int[] b = bonds[j];
                double r2 = sites.pos[b[0]].Mv1Squared(sites.pos[b[1]]);
                System.out.println(j + " " + b[0] + " " + b[1] + " " + Math.sqrt(r2));
            }
            for (int j = 0; j < angles.length; j++) {
                int[] b = angles[j];
                Vector v01 = new Vector3D();
                Vector v21 = new Vector3D();
                v01.Ev1Mv2(sites.pos[b[0]], sites.pos[b[1]]);
                v21.Ev1Mv2(sites.pos[b[2]], sites.pos[b[1]]);
                double dot = v21.dot(v01);
                double cos = dot / Math.sqrt(v01.squared() * v21.squared());
                double theta = 180*Math.acos(cos)/Math.PI;
                System.out.println(j + " " + b[0] + " " + b[1] + " " + b[2] + " " + theta);
            }
        }
    }
}