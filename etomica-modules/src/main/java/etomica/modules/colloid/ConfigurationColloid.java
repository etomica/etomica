/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class ConfigurationColloid implements Configuration {

    public ConfigurationColloid(Space space, ISpecies species, ISpecies speciesColloid, IRandom random) {
        this.space = space;
        this.species = species;
        this.speciesColloid = speciesColloid;
        this.random = random;
    }
    
    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }
    
    public int getChainLength() {
        return chainLength;
    }
    
    public void setNGraft(int newNGraft) {
        if (newNGraft == nGraft) return;
        nGraft = newNGraft;
        graftVectors = new Vector[nGraft];
        for (int i=0; i<nGraft; i++) {
            graftVectors[i] = space.makeVector();
        }
        if (nGraft == 1) {
            graftVectors[0].setX(0,1);
        }
        else if (nGraft == 2) {
            graftVectors[0].setX(0,1);
            graftVectors[1].setX(0,-1);
        }
        else if (nGraft == 4) {
            double x = 1.0/Math.sqrt(3);
            graftVectors[0].E(new double[]{ x, x, x});
            graftVectors[1].E(new double[]{-x,-x, x});
            graftVectors[2].E(new double[]{-x, x,-x});
            graftVectors[3].E(new double[]{ x,-x,-x});
        }
        else if (nGraft == 6) {
            graftVectors[0].setX(0,1);
            graftVectors[1].setX(0,-1);
            graftVectors[2].setX(1,1);
            graftVectors[3].setX(1,-1);
            graftVectors[4].setX(2,1);
            graftVectors[5].setX(2,-1);
        }
        else if (nGraft == 8) {
            Vector v = space.makeVector();
            v.E(1);
            for (int i=0; i<8; i++) {
                graftVectors[i].E(1.0/Math.sqrt(3));
                graftVectors[i].TE(v);
                for (int j=2; j>-1; j--) {
                    if (v.getX(j) == 1) {
                        v.setX(j, -1);
                        break;
                    }
                    v.setX(j, 1);
                }
            }
        }
        else if (nGraft == 12) {
            double phi = (1+Math.sqrt(5))/2;
            double x1 = 1/Math.sqrt(phi*phi+1);
            double x2 = phi*x1;

            graftVectors[0].E(new double[]{  0, x1, x2});
            graftVectors[1].E(new double[]{  0, x1,-x2});
            graftVectors[2].E(new double[]{  0,-x1, x2});
            graftVectors[3].E(new double[]{  0,-x1,-x2});
            graftVectors[4].E(new double[]{ x1, x2,  0});
            graftVectors[5].E(new double[]{ x1,-x2,  0});
            graftVectors[6].E(new double[]{-x1, x2,  0});
            graftVectors[7].E(new double[]{-x1,-x2,  0});
            graftVectors[8].E(new double[]{ x2,  0, x1});
            graftVectors[9].E(new double[]{ x2,  0,-x1});
            graftVectors[10].E(new double[]{-x2,  0, x1});
            graftVectors[11].E(new double[]{-x2,  0,-x1});
        }
        else if (nGraft == 20) {
            double phi = (1+Math.sqrt(5))/2;
            double x1 = 1/Math.sqrt(3);
            double x2 = x1*phi;
            double x3 = x1/phi;
            graftVectors[0].E(new double[]{ x1, x1, x1});
            graftVectors[1].E(new double[]{ x1, x1,-x1});
            graftVectors[2].E(new double[]{ x1,-x1, x1});
            graftVectors[3].E(new double[]{ x1,-x1,-x1});
            graftVectors[4].E(new double[]{-x1, x1, x1});
            graftVectors[5].E(new double[]{-x1, x1,-x1});
            graftVectors[6].E(new double[]{-x1,-x1, x1});
            graftVectors[7].E(new double[]{-x1,-x1,-x1});
            graftVectors[8].E(new double[]{  0, x3, x2});
            graftVectors[9].E(new double[]{  0, x3,-x2});
            graftVectors[10].E(new double[]{  0,-x3, x2});
            graftVectors[11].E(new double[]{  0,-x3,-x2});
            graftVectors[12].E(new double[]{ x3, x2,  0});
            graftVectors[13].E(new double[]{ x3,-x2,  0});
            graftVectors[14].E(new double[]{-x3, x2,  0});
            graftVectors[15].E(new double[]{-x3,-x2,  0});
            graftVectors[16].E(new double[]{ x2,  0, x3});
            graftVectors[17].E(new double[]{ x2,  0,-x3});
            graftVectors[18].E(new double[]{-x2,  0, x3});
            graftVectors[19].E(new double[]{-x2,  0,-x3});
        }
        else {
            throw new RuntimeException("don't know how to do nGraft = "+nGraft);
        }
    }
    
    public int getNGraft() {
        return nGraft;
    }
    
    public void setSigmaColloid(double newSigmaColloid) {
        sigmaColloid = newSigmaColloid;
    }
    
    public void setSigmaMonomer(double newSigmaMonomer) {
        sigmaMonomer = newSigmaMonomer;
    }
    
    public double getSigmaColloid() {
        return sigmaColloid;
    }
    
    public double getSigmaMonomer() {
        return sigmaMonomer;
    }
    
    public void setMonomerMonomerBondManager(AtomLeafAgentManager newMonomerMonomerBondManager) {
        monomerMonomerBondManager = newMonomerMonomerBondManager;
    }
    
    public void setColloidMonomerBondManager(AtomLeafAgentManager newColloidMonomerBondManager) {
        colloidMonomerBondManager = newColloidMonomerBondManager;
    }
    
    public void initializeCoordinates(Box box) {
        box.setNMolecules(species, 0);
        box.setNMolecules(species, nGraft*chainLength);
        box.setNMolecules(speciesColloid, 1);

        IMoleculeList colloidList = box.getMoleculeList(speciesColloid);
        IAtom colloidAtom = colloidList.getMolecule(0).getChildList().getAtom(0);
        ((AtomArrayList)colloidMonomerBondManager.getAgent(colloidAtom)).clear();
        Vector colloidPos = colloidAtom.getPosition();
        colloidPos.E(0);

        if (chainLength*nGraft == 0) {
            return;
        }
        
        IMoleculeList monomerList = box.getMoleculeList(species);
        IAtom previousAtom = null;
        int iMonomer = 0;
        Vector dr = space.makeVector();
        Vector temp = space.makeVector();
        Vector lastPos = space.makeVector();

        for (int j=0; j<nGraft; j++) {
            lastPos.E(colloidPos);
            lastPos.PEa1Tv1(0.5*sigmaColloid-0.5*sigmaMonomer, graftVectors[j]);
            for (int k=0; k<chainLength; k++) {
                IAtom atom = monomerList.getMolecule(iMonomer).getChildList().getAtom(0);
                if (k==0) {
                    ((AtomArrayList)colloidMonomerBondManager.getAgent(colloidAtom)).add(atom);
                    ((AtomArrayList)colloidMonomerBondManager.getAgent(atom)).add(colloidAtom);
                }
                else {
                    ((AtomArrayList)monomerMonomerBondManager.getAgent(previousAtom)).add(atom);
                    ((AtomArrayList)monomerMonomerBondManager.getAgent(atom)).add(previousAtom);
                }

                dr.E(graftVectors[j]);

                if (k>0) {
                    double tempSq = 0;
                    do {
                        // first get a random unit vector
                        temp.setRandomSphere(random);
                        // find the component of the unit vector perpendicular to our direction
                        temp.PEa1Tv1(-temp.dot(dr), dr);
                        // if the random unit vector was nearly parallel (or anti-parallel)
                        // to direction then the calculations will not be particularly
                        // precise, so try again
                        tempSq = temp.squared();
                    } while (tempSq < 0.001);
                    temp.TE(1/Math.sqrt(tempSq));
                    dr.TE(Math.cos(Math.PI/4));
                    dr.PEa1Tv1(Math.sin(Math.PI/4), temp);
                }

                Vector p = atom.getPosition();
                p.E(lastPos);
                p.PEa1Tv1(0.99*sigmaMonomer, dr);
                iMonomer++;
                previousAtom = atom;
                lastPos.E(p);
            }
        }
    }

    protected final Space space;
    protected ISpecies species, speciesColloid;
    protected int chainLength, nGraft;
    protected Vector[] graftVectors;
    protected double sigmaColloid, sigmaMonomer;
    protected AtomLeafAgentManager monomerMonomerBondManager, colloidMonomerBondManager;
    protected final IRandom random;
}
