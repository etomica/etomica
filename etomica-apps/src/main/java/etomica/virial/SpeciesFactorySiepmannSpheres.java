/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.config.ConformationChainZigZag;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

/**
 * SpeciesFactory that makes Siepmann's alkane model.
 */
public class SpeciesFactorySiepmannSpheres implements SpeciesFactory, java.io.Serializable {

    public SpeciesFactorySiepmannSpheres(Space space, int nA) {
        this(space, nA, nominalBondL, nominalBondTheta);
    }
    
    public SpeciesFactorySiepmannSpheres(Space space, int nA, double bondL, double bondTheta) {
        this.nA = nA;
        this.bondL = bondL;
        this.bondTheta = bondTheta;
        this.space = space;
        init();
    }
    
    public void setBondL(double newBondL) {
        bondL = newBondL;
        init();
    }
    
    public double getBondL() {
        return bondL;
    }
    
    public void setBondTheta(double newBondTheta) {
        bondTheta = newBondTheta;
        init();
    }
    
    public double getBondTheta() {
        return bondTheta;
    }
    
    public void init() {
        Vector vector1 = space.makeVector();
        vector1.setX(0, bondL);
        Vector vector2 = space.makeVector();
        vector2.setX(0, -bondL*Math.cos(bondTheta));
        vector2.setX(1, bondL*Math.sin(bondTheta));
        conformation = new ConformationChainZigZag(space, vector1, vector2);
    }
    
    public ISpecies makeSpecies(Space _space) {
        SpeciesAlkane species = new SpeciesAlkane(_space, nA);
        species.setConformation(conformation);
        return species;
    }
    
    private static final long serialVersionUID = 1L;
    protected static final double nominalBondL = 1.54;
    protected static final double nominalBondTheta = Math.PI*114/180;
    protected final Space space;
    protected double bondL;
    protected double bondTheta;
    private final int nA;
    private ConformationChainZigZag conformation;
}
