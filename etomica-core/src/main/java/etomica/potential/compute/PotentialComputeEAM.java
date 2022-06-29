/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.NeighborIteratorCell;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialEmbedding;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;
import etomica.util.collections.IntSet;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

/**
 * PotentialCompute that handles EAM systems consisting of pair and embedding
 * potentials.  This class can handle purely pairwise systems, although less
 * efficiently than PotentialComputePair.
 */
public class PotentialComputeEAM implements PotentialCompute {

    private final NeighborIterator neighborIterator;
    protected final IPotentialEmbedding[] embeddingPotentials;
    protected final IPotential2[][] pairPotentials;
    protected final IPotential2[] rhoPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    private final NeighborIteratorCell.SuperNbrConsumer nbrConsumer1, nbrConsumerEmbedding;
    protected double[] rhoSum, idf, id2f;
    protected double[] uAtom;
    protected final DoubleArrayList rdrho, drhoSum;
    protected final HashMap<IntSet,PairData> pairDataHashMap;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged, rhoAtomsChanged;
    protected final Vector zero;
    protected double virialTot = Double.NaN, virialTot2 = Double.NaN, energyTot = Double.NaN;
    protected double[] virialXYZ = new double[6];
    protected double[] virial2XYZ = new double[9];
    protected Vector[] forces , dFdeV;
    protected Vector[][] dFde;
    protected final Space space;
    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;
    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;
    public Tensor unity, Hij, Hik, Hkj;

    public PotentialComputeEAM(SpeciesManager sm, Box box, NeighborManager neighborManager) {
        space = box.getSpace();
        int typeCount = sm.getAtomTypeCount();
        pairPotentials = new IPotential2[typeCount][typeCount];
        rhoPotentials = new IPotential2[typeCount];
        embeddingPotentials = new IPotentialEmbedding[typeCount];
        this.neighborManager = neighborManager;
        this.neighborManager.setPairPotentials(pairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        rhoAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        rdrho = new DoubleArrayList(16);
        drhoSum = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        Hij = space.makeTensor();
        Hik = space.makeTensor();
        Hkj = space.makeTensor();
        unity = space.makeTensor();
        unity.setComponent(0, 0, 1.0);
        unity.setComponent(1, 1, 1.0);
        unity.setComponent(2, 2, 1.0);

        forces = new Vector[0];
        dFdeV = new Vector[0];
        dFde  = new Vector[9][0];

        pairDataHashMap = new HashMap<>();

        this.atomCountByType = new int[typeCount];
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]++;
                }
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });

        // for new atom energy, includes embedding contributions
        this.nbrConsumer1 = new NeighborIteratorCell.SuperNbrConsumer() {
            @Override
            public double accept(IAtom iAtom, IAtom jAtom, Vector rij, int n) {
                int i = iAtom.getLeafIndex();
                int j = jAtom.getLeafIndex();
                int iType = iAtom.getType().getIndex();
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = pairPotentials[iType][jType];
                double r2 = rij.squared();
                double uij = 0;
                if (pij != null) {
                    uij += pij.u(r2);
                    uAtomsChanged.add(jAtom.getLeafIndex());
                    if (duAtomMulti) {
                        duAtom.plusEquals(i, 0.5 * uij);
                        duAtom.plusEquals(j, 0.5 * uij);
                    } else {
                        duAtom.plusEquals(0, 0.5 * uij);
                        duAtom.add(0.5 * uij);
                    }
                }
                IPotential2 irp = rhoPotentials[iType];
                IPotential2 jrp = rhoPotentials[jType];
                if (jrp != null) {
                    double jrc = jrp.getRange();
                    if (r2 <= jrc * jrc) {
                        double rho = jrp.u(r2);
                        drhoSum.plusEquals(i, rho);
                        // if i doesn't contribute to j's density, then bail
                        if (iType != jType) {
                            if (irp == null) return uij;
                            double irc = irp.getRange();
                            if (r2 > irc * irc) return uij;
                            rho = irp.u(r2);
                        }
                        rhoAtomsChanged.add(j);
                        // we need the energy of the configuration with the atom in the new spot
                        // minus the energy with no atom
                        uij += embeddingPotentials[jType].u(rhoSum[j] + drhoSum.getDouble(j) + rho);
                        uij -= embeddingPotentials[jType].u(rhoSum[j] + drhoSum.getDouble(j));
                        // drhoSum will hold new-old
                        // processAtom(+1) handles this and we should be done
                        // we'll be called again and recompute rho for the old config
                        // for processAtom(-1), we ignore drhoSum since it was already handled
                        drhoSum.plusEquals(j, rho);
                    }
                }
                return uij;
            }
        };

        // for old embedding energy
        this.nbrConsumerEmbedding = new NeighborIteratorCell.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom jAtom, Vector rij, int n) {
                double r2 = rij.squared();
                int jType = jAtom.getType().getIndex();
                double jrc = rhoPotentials[jType].getRange();
                if (r2 > jrc * jrc) return 0;

                int j = jAtom.getLeafIndex();
                rhoAtomsChanged.add(j);
                double rho = rhoPotentials[jType].u(r2);
                drhoSum.replace(j, -rho);
                // we need to compute the old energy of the system with and without iAtom
                return embeddingPotentials[jType].u(rhoSum[j])
                        - embeddingPotentials[jType].u(rhoSum[j] - rho);
            }
        };
    }

    @Override
    public boolean needForcesForVirial() {
        return true;
    }

    @Override
    public void init() {
        this.neighborManager.init();
    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    public Vector[] getdFdeV()  { return dFdeV; }

    public Vector[][] getdFde()  { return dFde; }

    @Override
    public double getLastVirial() {
        return virialTot;
    }

    public double getLastVirial2() {
        return virialTot2;
    }

    public double[] getLastVirialXYZ() {
        return virialXYZ;
    }

    public double[] getLastVirial2XYZ() {
        return virial2XYZ;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, IPotential2 p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;

        updateMaxRange();
    }

    public void setRhoPotential(AtomType atomType, IPotential2 pRho) {
        rhoPotentials[atomType.getIndex()] = pRho;

        updateMaxRange();
    }

    private void updateMaxRange() {

        double maxRange1 = Arrays.stream(pairPotentials).flatMap(Arrays::stream)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential2::getRange)
                .max().orElse(0);

        double maxRange2 = Arrays.stream(rhoPotentials)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential2::getRange)
                .max().orElse(0);

        this.neighborManager.setPotentialRange(Math.max(maxRange1, maxRange2));
    }

    public void setEmbeddingPotential(AtomType atomType, IPotentialEmbedding p) {
        embeddingPotentials[atomType.getIndex()] = p;
    }

    @Override
    public void updateAtom(IAtom atom) {
        this.neighborManager.updateAtom(atom);
    }

    protected final void zeroArrays(boolean doForces, boolean wantsHessian) {
        virialTot = 0;
        virialTot2 = 0;
        for (int i=0; i<6; i++) virialXYZ[i] = 0;
        for (int i=0; i<9; i++) virial2XYZ[i] = 0;

        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
        }

        if (wantsHessian && numAtoms > dFdeV.length) {
            int oldLength = dFdeV.length;
            dFdeV = Arrays.copyOf(dFdeV, numAtoms);
            dFde  = new Vector[9][numAtoms];


            for (int i = oldLength; i < numAtoms; i++) dFdeV[i] = box.getSpace().makeVector();
            for (int i=0; i < 9; i++){
                for (int j = oldLength; j < numAtoms; j++){
                    dFde[i][j] = box.getSpace().makeVector();
                }
            }
        }

        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
            rhoSum = new double[numAtoms];
            idf = new double[numAtoms];
            id2f = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            rhoSum[i] = 0;
            if (doForces) forces[i].E(0);
            if (wantsHessian) {
                dFdeV[i].E(0);
                for (int j=0; j<9; j++){
                    dFde[j][i].E(0);
                }
            }
        }
    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        final boolean wantsHessian = pc != null && pc.wantsHessian();
        zeroArrays(doForces, wantsHessian);
        rdrho.clear();
        pairDataHashMap.clear();

        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.size();

        double[] uTot = {0};

        for (int i = 0; i < numAtoms; i++) { //i
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotential2[] ip = pairPotentials[iType];
            IPotential2 irp = rhoPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(finalI, (jAtom, rij, n) -> { //j
                int j = jAtom.getLeafIndex();
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                double[] u012 = new double[3];
                double r2 = rij.squared();

                // Pairwise contribution
                if (pij != null) {
                    pij.u012add(r2, u012);
                    double uij = u012[0];
                    if (uij != 0) { //or r2>rc
                        uAtom[finalI] += 0.5 * uij;
                        uAtom[j] += 0.5 * uij;
                        double duij = u012[1];

                        if (doForces) {
                            virialTot += duij; //Fij.rij
                            virialXYZ[0] += duij/r2*rij.getX(0)*rij.getX(0); //xx
                            virialXYZ[1] += duij/r2*rij.getX(1)*rij.getX(1); //yy
                            virialXYZ[2] += duij/r2*rij.getX(2)*rij.getX(2); //zz
                            virialXYZ[3] += duij/r2*rij.getX(1)*rij.getX(2); //yz
                            virialXYZ[4] += duij/r2*rij.getX(0)*rij.getX(2); //xz
                            virialXYZ[5] += duij/r2*rij.getX(0)*rij.getX(1); //xy

                            forces[finalI].PEa1Tv1(duij/r2, rij);
                            forces[j].PEa1Tv1(-duij/r2, rij);
                        }
                        uTot[0] += uij;

                        if (wantsHessian) {
                            pc.pairCompute(finalI, j, rij, u012);
                            double d2uij = u012[2];
                            virialTot2 += 1.0/9.0*d2uij;
                            double fac = (d2uij-duij)/r2/r2;
                            virial2XYZ[0] += fac*rij.getX(0)*rij.getX(0)*rij.getX(0)*rij.getX(0);//C11
                            virial2XYZ[1] += fac*rij.getX(1)*rij.getX(1)*rij.getX(1)*rij.getX(1);//C22
                            virial2XYZ[2] += fac*rij.getX(2)*rij.getX(2)*rij.getX(2)*rij.getX(2);//C33
                            virial2XYZ[3] += fac*rij.getX(0)*rij.getX(0)*rij.getX(1)*rij.getX(1);//C12
                            virial2XYZ[4] += fac*rij.getX(0)*rij.getX(0)*rij.getX(2)*rij.getX(2);//C13
                            virial2XYZ[5] += fac*rij.getX(1)*rij.getX(1)*rij.getX(2)*rij.getX(2);//C23
                            virial2XYZ[6] += fac*rij.getX(1)*rij.getX(2)*rij.getX(1)*rij.getX(2);//C44
                            virial2XYZ[7] += fac*rij.getX(0)*rij.getX(2)*rij.getX(0)*rij.getX(2);//C55
                            virial2XYZ[8] += fac*rij.getX(0)*rij.getX(1)*rij.getX(0)*rij.getX(1);//C66

                            dFdeV[finalI].PEa1Tv1(d2uij/r2/3.0, rij);
                            dFdeV[j].PEa1Tv1(-d2uij/r2/3.0, rij);

                            Vector tmpVec = space.makeVector();
                            int indCount = 0;
                            //00 , 01 , 02 , 11 , 12 , 22
                            //0    1    2    3    4    5
                            for (int l=0; l<3; l++){
                                for (int ll=l; ll<3; ll++){
                                    tmpVec.E(0);
                                    if (l==ll) {
                                        tmpVec.setX(l, duij/r2*rij.getX(l));
                                    } else {
                                        tmpVec.setX(l,  0.5*duij/r2*rij.getX(ll));
                                        tmpVec.setX(ll, 0.5*duij/r2*rij.getX(l));
                                    }
                                    tmpVec.PEa1Tv1(fac*rij.getX(l)*rij.getX(ll), rij);
                                    dFde[indCount][finalI].PE(tmpVec);
                                    dFde[indCount][j].ME(tmpVec);
                                    indCount++;
                                }
                            }

                        }

                    }
                    u012[0] = u012[1] = u012[2] = 0;
                }

                // Multi-body contribution
                if (irp != null && embeddingPotentials[jType] != null) {
                    // i contributes to j density: rho_i
                    double irc = irp.getRange();
                    if (r2 <= irc * irc) {
                        u012[0] = u012[1] = u012[2] = 0;
                        irp.u012add(r2, u012);
                        rdrho.add(u012[1]);
                        rhoSum[j] += u012[0];
                        if (wantsHessian) {
                            PairData pd = new PairData(u012[1], u012[2]);
                            if (iType == jType) {
                                int i1 = Math.min(finalI, j);
                                pairDataHashMap.put(new IntSet(i1, finalI + j - i1), pd);
                            } else {
                                pairDataHashMap.put(new IntSet(finalI, j), pd);
                            }
                        }
                    }

                    if (iType != jType) {
                        IPotential2 jrp = rhoPotentials[jType];
                        if (jrp != null && embeddingPotentials[iType] != null) {
                            double jrc = jrp.getRange();
                            if (r2 <= jrc * jrc) {
                                // j contributes to i density: rho_j
                                u012[0] = u012[1] = u012[2] = 0;
                                jrp.u012add(r2, u012);
                                rdrho.add(u012[1]);
                                rhoSum[finalI] += u012[0];
                                if (wantsHessian) {
                                    PairData pd = new PairData(u012[1], u012[2]);
                                    pairDataHashMap.put(new IntSet(j, finalI), pd);
                                }
                            }
                        }
                    } else {
                        rhoSum[finalI] += u012[0];
                    }
                } else if (iType != jType) {
                    IPotential2 jrp = rhoPotentials[jType];
                    if (jrp != null && embeddingPotentials[iType] != null) {
                        double jrc = jrp.getRange();
                        if (r2 <= jrc * jrc) {
                            // j contributes to i density
                            u012[0] = u012[1] = u012[2] = 0;
                            jrp.u012add(r2, u012);
                            rdrho.add(u012[1]);
                            rhoSum[finalI] += u012[0];
                            if (wantsHessian) {
                                PairData pd = new PairData(u012[1], u012[2]);
                                pairDataHashMap.put(new IntSet(j, finalI), pd);
                            }
                        }
                    }
                }//multi-body


            });//j loop
        }// i loop

        //uTot (pair and multibody) and dF/dtheta
        for (int i = 0; i < numAtoms; i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            if (embeddingPotentials[iType] == null) continue;
            double[] u = {0}, du = {0}, d2u = {0};
//            System.out.println(ElectronVolt.UNIT.fromSim(Math.sqrt(rhoSum[i]))*ElectronVolt.UNIT.fromSim(Math.sqrt(rhoSum[i])));
            embeddingPotentials[iType].udud2u(rhoSum[i], u, du, d2u);
//            System.out.println(ElectronVolt.UNIT.fromSim(u[0]));

            uTot[0] += u[0];
            if (doForces) {
                idf[i] = du[0]; //idf = dF_i/dtheta_i
                id2f[i] = d2u[0]; //idf = d2F_i/dtheta_i^2
            }
        }

        int[] rdrhoIdx = {0};

        if (doForces) {
            for (int i = 0; i < numAtoms; i++) {
                IAtom iAtom = atoms.get(i);
                int iType = iAtom.getType().getIndex();

                IPotential2 irp = rhoPotentials[iType];
                double irc = irp.getRange();
                int finalI = i;
                neighborIterator.iterUpNeighbors(i, (jAtom, rij, n) -> {
                    int j = jAtom.getLeafIndex();
                    int jType = jAtom.getType().getIndex();
                    if (rhoPotentials[jType] == null) return;
                    double r2 = rij.squared();

                    if (r2 <= irc*irc) {
                        double fac;
                        if (iType == jType) {
                            fac = (idf[finalI] + idf[j]) * rdrho.getDouble(rdrhoIdx[0]);
                            rdrhoIdx[0]++;
                        } else {
                            fac = idf[finalI] * rdrho.getDouble(rdrhoIdx[0]);
                            rdrhoIdx[0]++;
                            double jrc = rhoPotentials[jType].getRange();
                            if (r2 <= jrc * jrc) {
                                fac += idf[j] * rdrho.getDouble(rdrhoIdx[0]);
                                rdrhoIdx[0]++;
                            }
                        }
                        virialTot += fac; //Fij.rij
                        virialXYZ[0] += fac/r2*rij.getX(0)*rij.getX(0);//xx
                        virialXYZ[1] += fac/r2*rij.getX(1)*rij.getX(1);//yy
                        virialXYZ[2] += fac/r2*rij.getX(2)*rij.getX(2);//zz
                        virialXYZ[3] += fac/r2*rij.getX(1)*rij.getX(2);//yz
                        virialXYZ[4] += fac/r2*rij.getX(0)*rij.getX(2);//xz
                        virialXYZ[5] += fac/r2*rij.getX(0)*rij.getX(1);//xy

                        fac /= r2;
                        rij.TE(fac);
                        forces[finalI].PE(rij);
                        forces[j].ME(rij);
                    }
                });
            }
        }
        rdrhoIdx[0] = 0;

        if (wantsHessian) {
            virialTot2 -= 2.0/9.0*virialTot; //Fij.rij (pair + embedding)
            Vector tmpV = space.makeVector();
            for (int k = 0; k < numAtoms; k++) {
                IAtom kAtom = atoms.get(k);
                int kType = kAtom.getType().getIndex();
                IPotential2 krp = rhoPotentials[kType];
                double krc = krp.getRange();
                int finalK = k;
                neighborIterator.iterAllNeighbors(finalK, (iAtom, rki, n) -> {
                    double rki2 = rki.squared();
                    int iType = iAtom.getType().getIndex();
                    if (rhoPotentials[iType] == null || rki2 > krc*krc) return;
                    int i = iAtom.getLeafIndex();

                    PairData pdik, pdki;
                    if (iType == kType) {
                        int ik1 = Math.min(finalK, i);
                        pdik = pairDataHashMap.get(new
                                IntSet(ik1, finalK+i-ik1));
                        pdki = pdik;
                    } else {
                        pdik = pairDataHashMap.get(new IntSet(i, finalK));
                        pdki = pairDataHashMap.get(new IntSet(finalK, i));
                    }
                    double rdrhoik  = pdik.rdrho;

                    // Direct Hik
                    double rdrho2ik = pdik.r2drho;
                    double rdrhoki  = pdki.rdrho;
                    double rdrho2ki = pdki.r2drho;

                    double fac;
                    if (i < finalK) {
                        fac = (idf[i]*(rdrhoki-rdrho2ki)+idf[finalK]*(rdrhoik-rdrho2ik)-id2f[i]*rdrhoki*rdrhoki-id2f[finalK]*rdrhoik*rdrhoik)/rki2/rki2;
                        Hik.Ev1v2(rki, rki);
                        Hik.TE(fac);
                        fac = -(idf[i]*rdrhoki + idf[finalK]*rdrhoik)/rki2;
                        Hik.PEa1Tt1(fac, unity);
                        pc.pairComputeHessian(i,finalK, Hik);
                    }

                    // dF/deV
                    // Hik = d2uk/drki drkj
                    Hik.Ev1v2(rki, rki);
                    fac = (idf[finalK]*(rdrho2ik-rdrhoik)+id2f[finalK]*rdrhoik*rdrhoik)/rki2/rki2;
                    Hik.TE(fac);
                    double fac1 = idf[finalK]*rdrhoik/rki2;
                    Hik.PEa1Tt1(fac1, unity);
                    tmpV.E(rki);
                    Hik.transform(tmpV);
                    tmpV.TE(1.0/3.0);
                    dFdeV[finalK].PE(tmpV);
                    dFdeV[i].ME(tmpV);

                    //dF/dex
                    Vector tmpVec = space.makeVector();
                    int indCount = 0;
                    //00 , 01 , 02 , 11 , 12 , 22
                    //0    1    2    3    4    5
                    for (int l=0; l<3; l++){
                        for (int ll=l; ll<3; ll++){
                            tmpVec.E(0);
                            if(l==ll) {
                                tmpVec.setX(ll,fac1*rki.getX(l));
                            } else {
                                tmpVec.setX(l, 0.5*fac1*rki.getX(ll));
                                tmpVec.setX(ll,0.5*fac1*rki.getX(l));
                            }
                            tmpVec.PEa1Tv1(fac*rki.getX(l)*rki.getX(ll), rki);
                            dFde[indCount][finalK].PE(tmpVec);
                            dFde[indCount][i].ME(tmpVec);
                            indCount++;
                        }
                    }



                    //2nd virial
                    if(i < finalK) {
                        virialTot2 += 2.0 / 3.0 * tmpV.dot(rki);
                        fac  = (idf[finalK] * (rdrho2ik - rdrhoik) + id2f[finalK] * rdrhoik * rdrhoik
                              + idf[i] * (rdrho2ki - rdrhoki)      + id2f[i] * rdrhoki * rdrhoki)/rki2/rki2;
                        virial2XYZ[0] += fac * rki.getX(0) * rki.getX(0) * rki.getX(0) * rki.getX(0);
                        virial2XYZ[1] += fac * rki.getX(1) * rki.getX(1) * rki.getX(1) * rki.getX(1);
                        virial2XYZ[2] += fac * rki.getX(2) * rki.getX(2) * rki.getX(2) * rki.getX(2);
                        virial2XYZ[3] += fac * rki.getX(0) * rki.getX(0) * rki.getX(1) * rki.getX(1);
                        virial2XYZ[4] += fac * rki.getX(0) * rki.getX(0) * rki.getX(2) * rki.getX(2);
                        virial2XYZ[5] += fac * rki.getX(1) * rki.getX(1) * rki.getX(2) * rki.getX(2);
                        virial2XYZ[6] += fac * rki.getX(1) * rki.getX(2) * rki.getX(1) * rki.getX(2);
                        virial2XYZ[7] += fac * rki.getX(0) * rki.getX(2) * rki.getX(0) * rki.getX(2);
                        virial2XYZ[8] += fac * rki.getX(0) * rki.getX(1) * rki.getX(0) * rki.getX(1);
                    }
                    neighborIterator.iterAllNeighbors(finalK, (jAtom, rkj, m) -> {
                        int j = jAtom.getLeafIndex();
                        int jType = jAtom.getType().getIndex();
                        double rkj2 = rkj.squared();
                        if (j <= i || rhoPotentials[jType] == null || rkj2 > krc*krc) return;
                        PairData pdjk;
                        if (jType == kType) {
                            int jk1 = Math.min(j, finalK);
                            pdjk = pairDataHashMap.get(new IntSet(jk1, j+finalK-jk1));
                        } else {
                            pdjk = pairDataHashMap.get(new IntSet(j, finalK));
                        }

                        double rdrhojk  = pdjk.rdrho;

                        //Indirect Hij, Hik, Hkj
                        Hij.Ev1v2(rki, rkj);
                        Hij.TE(id2f[finalK]*rdrhoik*rdrhojk/rki2/rkj2);
                        pc.pairComputeHessian(i,j,Hij);

                        Hik.E(Hij);
                        Hik.TE(-1);
                        pc.pairComputeHessian(i,finalK,Hik);

                        Hkj.E(Hij);
                        Hkj.TE(-1);
                        pc.pairComputeHessian(finalK,j,Hkj);


                        //dF/deV
                        //ki
                        Hij.Ev1v2(rki, rkj);
                        double fac2 = id2f[finalK]*rdrhoik*rdrhojk/rki2/rkj2;
                        Hij.TE(fac2);
                        tmpV.E(rkj);
                        Hij.transform(tmpV);
                        tmpV.TE(1.0/3.0);
                        dFdeV[finalK].PE(tmpV);
                        dFdeV[i].ME(tmpV);
                        //kj
                        tmpV.E(rki);
                        Hij.transpose();
                        Hij.transform(tmpV);
                        tmpV.TE(1.0/3.0);
                        dFdeV[finalK].PE(tmpV);
                        dFdeV[j].ME(tmpV);

                        //dF/dex
                        int indCount2 = 0;
                        //00 , 01 , 02 , 11 , 12 , 22
                        //0    1    2    3    4    5
                        for (int l=0; l<3; l++){
                            for (int ll=l; ll<3; ll++){
                                //ki
                                dFde[indCount2][finalK].PEa1Tv1(fac2*rkj.getX(l)*rkj.getX(ll), rki);
                                dFde[indCount2][i].PEa1Tv1(-fac2*rkj.getX(l)*rkj.getX(ll), rki);
                                //kj
                                dFde[indCount2][finalK].PEa1Tv1(fac2*rki.getX(l)*rki.getX(ll), rkj);
                                dFde[indCount2][j].PEa1Tv1(-fac2*rki.getX(l)*rki.getX(ll), rkj);
                                indCount2++;
                            }
                        }




                        //2nd virial
                        virialTot2 += 2.0/3.0*tmpV.dot(rkj);
                        virial2XYZ[0] += 2*fac2*rki.getX(0)*rki.getX(0)*rkj.getX(0)*rkj.getX(0);
                        virial2XYZ[1] += 2*fac2*rki.getX(1)*rki.getX(1)*rkj.getX(1)*rkj.getX(1);
                        virial2XYZ[2] += 2*fac2*rki.getX(2)*rki.getX(2)*rkj.getX(2)*rkj.getX(2);
                        virial2XYZ[3] += 2*fac2*rki.getX(0)*rki.getX(0)*rkj.getX(1)*rkj.getX(1);
                        virial2XYZ[4] += 2*fac2*rki.getX(0)*rki.getX(0)*rkj.getX(2)*rkj.getX(2);
                        virial2XYZ[5] += 2*fac2*rki.getX(1)*rki.getX(1)*rkj.getX(2)*rkj.getX(2);
                        virial2XYZ[6] += 2*fac2*rki.getX(1)*rki.getX(2)*rkj.getX(1)*rkj.getX(2);
                        virial2XYZ[7] += 2*fac2*rki.getX(0)*rki.getX(2)*rkj.getX(0)*rkj.getX(2);
                        virial2XYZ[8] += 2*fac2*rki.getX(0)*rki.getX(1)*rkj.getX(0)*rkj.getX(1);
                    });//j
                });//i
            }//k
        }//if (wantsHessian)

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot[0] += uCorrection[0];
        virialTot += duCorrection[0];
        energyTot = uTot[0];
        return uTot[0];
    }//computeAll

    public double oldEmbeddingEnergy(IAtom iAtom) {
        // just compute all the embedding energies
        int i = iAtom.getLeafIndex();
        rhoAtomsChanged.clear();
        rhoAtomsChanged.add(i);
        int iType = iAtom.getType().getIndex();
        double u = embeddingPotentials[iType].u(rhoSum[i]);
        u += this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumerEmbedding);
        return u;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2 + oldEmbeddingEnergy(iAtom);
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += this.computeOneTruncationCorrection(atom.getLeafIndex());
        }
        return u;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        this.duAtomMulti = false;
        int i = iAtom.getLeafIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        uAtomsChanged.add(i);
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        duAtom.add(0);
        drhoSum.clear();
        drhoSum.ensureCapacity(atoms.size());
        if (rhoAtomsChanged.size() == 0) {
            // if we called oldEnergy for this atom, then it will already be in the list
            rhoAtomsChanged.add(i);
        }
        double u = this.computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(i);
        return u;
    }

    protected double computeOneInternal(IAtom iAtom) {
        return this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumer1);
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        throw new RuntimeException("EAM can't handle many atoms");
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        throw new RuntimeException("EAM can't handle many atoms");
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        duAtomMulti = true;
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (IAtom atom : molecule.getChildList()) {
            uAtomsChanged.add(atom.getLeafIndex());

            u += computeOneInternal(atom);
            u += computeOneTruncationCorrection(atom.getLeafIndex());
        }

        return u;
    }

    @Override
    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            if (duAtomMulti) {
                uAtom[jj] += fac * duAtom.getDouble(jj);
                duAtom.replace(jj, 0);
            } else {
                uAtom[jj] += fac * duAtom.getDouble(j);
            }
        }

        for (int j = 0; j < rhoAtomsChanged.size(); j++) {
            int jj = rhoAtomsChanged.getInt(j);
            // by the time we get to processAtom(+1), we have
            // the difference.  just use that
            if (fac == 1) {
                rhoSum[jj] += drhoSum.getDouble(jj);
            }
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this.neighborManager.makeIntegratorListener();
    }

    public void computeAllTruncationCorrection(double[] uCorrection, double[] duCorrection) {
        if (!doAllTruncationCorrection) {
            return;
        }
        double uCor = 0;
        double duCor = 0;
        for (int i = 0; i < atomCountByType.length; i++) {
            for (int j = i; j < atomCountByType.length; j++) {
                IPotential2 p = pairPotentials[i][j];
                int numPairs;
                if (j == i) {
                    numPairs = atomCountByType[i] * (atomCountByType[j] - 1) / 2;
                } else {
                    numPairs = atomCountByType[i] * atomCountByType[j];
                }
                double pairDensity = numPairs / box.getBoundary().volume();

                double[] u = new double[1];
                double[] du = new double[1];
                p.u01TruncationCorrection(space, u, du);
                uCor += pairDensity * u[0];
                duCor += pairDensity * du[0];

            }
        }
        uCorrection[0] = uCor;
        duCorrection[0] = duCor;
    }

    public double computeOneTruncationCorrection(int iAtom) {
        if (!doOneTruncationCorrection) {
            return 0;
        }
        int iType = box.getLeafList().get(iAtom).getType().getIndex();
        double uCorrection = 0;
        for (int j = 0; j < atomCountByType.length; j++) {
            IPotential2 p = pairPotentials[iType][j];
            double pairDensity;
            if (iType == j) {
                pairDensity = (atomCountByType[j] - 1) / box.getBoundary().volume();
            } else {
                pairDensity = atomCountByType[j] / box.getBoundary().volume();
            }
            double integral = p.integral(space, p.getRange());
            uCorrection += pairDensity * integral;
        }
        return uCorrection;
    }

    public static class PairData {
        public double rdrho, r2drho;
        public PairData(double rdrho, double r2drho) {
            this.rdrho = rdrho;
            this.r2drho = r2drho;
        }
    }
}