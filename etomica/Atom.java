package simulate;
import java.awt.*;

/**
*  Each instance of the class Atom holds the position and velocity of one
*  physical atom; all simulation kinetics and dynamics are performed by 
*  operating on these values.
*  
*  @author David Kofke
*  @author C. Daniel Barnes
*  @see Molecule
*/
public final class Atom implements Space.Occupant {

    /**
    * Constructs an atom with default values for mass, diameter, and color.
    * Default values are mass = 1.0 amu; diameter = 0.1 A; color = black.
    * Defaults for all coordinates and momenta are zero.
    *
    * @param parent       molecule in which atom resides
    * @param index        sequential index of atom as assigned by parent molecule
    */
    public Atom(Molecule parent, AtomType t, int index) {
        parentMolecule = parent;
        type = t;
        atomIndex = index;
        coordinate = parentMolecule.parentSpecies.parentSimulation.space.makeCoordinate(this);
        r = coordinate.position();
        p = coordinate.momentum();
        workVector = coordinate.makeVector();
        rLast = coordinate.makeVector();
        velocity = coordinate.makeVector();
        setStationary(false);
        useTypeColor();
    }
                    
    public void setIntegratorAgent(Integrator.Agent ia) {this.ia = ia;}
        
    public final Molecule parentMolecule() {return parentMolecule;}
        
    public final Space.Coordinate coordinate() {return coordinate;}
            
    public final int getSpeciesIndex() {return parentMolecule.getSpeciesIndex();}
    public final int getAtomIndex() {return atomIndex;}
        
    public final Color getColor() {return color;}
    public final void setColor(Color c) {this.color = c;}
    public final void useTypeColor() {this.color = type.color();}  //indicates that atom color is determined by its type
        
    public void setStationary(boolean b) {stationary = b;}
    public final boolean isStationary() {return stationary;}

    /**
    * Sets atom following this one in linked list, and sets this to be that
    * atom's previous atom in list
    * 
    * @param atom  the next atom in the list
    */
    public final void setNextAtom(Atom atom) {
        nextAtom = atom;
        if(atom!=null) {atom.previousAtom = this;}
    }
    public final void clearPreviousAtom() {previousAtom = null;}
    public final Atom nextAtom() {return nextAtom;}
    public final Atom previousAtom() {return previousAtom;} 
    
    public final Atom nextMoleculeFirstAtom() {return parentMolecule.lastAtom.nextAtom();}  //first atom on next molecule
    public final Atom previousMoleculeLastAtom() {return parentMolecule.firstAtom.previousAtom();}  //first atom on next molecule

    public final double mass() {return type.mass();}
    public final double rm() {return type.rm();}

    public final Phase parentPhase() {return parentMolecule.parentPhase;}

    public void draw(Graphics g, int[] origin, double scale) {type.draw(g, origin, scale, this);}

    public final void translateBy(Space.Vector u) {r.PE(u);}
    public final void translateBy(double d, Space.Vector u) {r.PEa1Tv1(d,u);}
    public final void translateTo(Space.Vector u) {r.E(u);}      
    public final void translateToRandom(simulate.Phase p) {translateTo(p.boundary().randomPosition());}
    public final void displaceBy(Space.Vector u) {rLast.E(r); translateBy(u);}
    public final void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,u);}
    public final void displaceTo(Space.Vector u) {rLast.E(r); translateTo(u);}  
    public final void displaceWithin(double d) {workVector.setRandomCube(); displaceBy(d,workVector);}
    public final void displaceToRandom(simulate.Phase p) {rLast.E(r); translateToRandom(p);}
    public final void replace() {r.E(rLast);}
    public final void inflate(double s) {r.TE(s);}

    public final void accelerateBy(Space.Vector u) {p.PE(u);}
    public final void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}

    public final double kineticEnergy() {return coordinate.kineticEnergy(type.mass());}
    public final void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
        double magnitude = Math.sqrt(type.mass()*temperature*(double)Simulation.D/Constants.KE2T);  //need to divide by sqrt(m) to get velocity
        p.setRandomSphere();
        p.TE(magnitude);
    }
    public final void scaleMomentum(double scale) {p.TE(scale);}

    public final Space.Vector position() {return r;}
    public final Space.Vector momentum() {return p;}
    public final double position(int i) {return r.component(i);}
    public final double momentum(int i) {return p.component(i);}
    public final Space.Vector velocity() {velocity.E(p); velocity.TE(type.rm()); return velocity;}  //returned vector is not thread-safe

    public Integrator.Agent ia;
        
    /**
    * Color of the atom when drawn on the screen
    * This color is set by the colorScheme object in the atom's species
    */
    Color color = Color.black;
        
    /**
    * Flag indicating whether atom is stationary or mobile.
    * Default is false (atom is mobile)
    */
    private boolean stationary;
       
    /**
    * Instance of molecule in which this atom resides.
    * Assigned in Atom constructor.
    * @see Molecule#makeAtoms
    */
//    final Molecule parentMolecule;
     Molecule parentMolecule;
        
    /**
    * Identifier of atom within molecule.
    * Assigned by parent molecule when invoking Atom constructor.
    * @see Molecule#makeAtoms
    */
//    final int atomIndex;
     int atomIndex;
    
    private Atom nextAtom, previousAtom;
        
    public  Space.Coordinate coordinate;
    public  Space.Vector r, p;  //position, momentum
        
    public  AtomType type;
//    public final Space.Coordinate coordinate;
//    public final Space.Vector r, p;  //position, momentum
        
//    public final AtomType type;
        
//    public final Space.Vector workVector, rLast, velocity;
    public  Space.Vector workVector, rLast, velocity;
        
    public interface Iterator {        
        public boolean hasNext();
        public Atom next();
        public void reset(Atom a);
        public void reset();
    
        public static final class Up implements Iterator {
            private Atom atom, nextAtom;
            private boolean hasNext;
            private Phase phase;
            public Up(Phase p) {phase = p; hasNext = false;}
            public Up(Phase p, Atom a) {phase = p; reset(a);}
            public boolean hasNext() {return hasNext;}
            public void reset(Atom a) {
                if(a == null) {hasNext = false; return;}
                atom = a;
                hasNext = true;
            }
            public void reset() {reset(phase.firstAtom());}
            public Atom next() {
                nextAtom = atom;
                atom = atom.nextAtom();
                if(atom == null) {hasNext = false;}
                return nextAtom;
            }
        } //end of Atom.Iterator.Up
        
        public static final class Down implements Iterator {
            private Atom atom, nextAtom, firstAtom;
            private boolean hasNext;
            private Phase phase;
            public Down(Phase p) {phase = p; hasNext = false;}
            public Down(Phase p, Atom a) {phase = p; reset(a);}
            public boolean hasNext() {return hasNext;}
            public void reset(Atom a) {
                if(a == null) {hasNext = false; return;}
                atom = a;
                hasNext = true;
            }
            public void reset() {reset(phase.firstAtom());}
            public Atom next() {
                nextAtom = atom;
                atom = atom.previousAtom();
                if(atom == null) {hasNext = false;}
                return nextAtom;
            }
        } //end of Atom.Iterator.Down
    }
}