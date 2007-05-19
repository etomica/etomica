package etomica.models.hexane;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.util.IRandom;


/**
 * MC move which performs a configuration bias monte carlo moved, then moves 
 * that atom to a new position so that its center of mass does not change.
 * 
 * @author cribbin
 *
 */
public class MCMoveCombinedCbmcTranslation extends MCMovePhase {

    protected MCMoveCBMC cbmcMove;
    protected double uNew, uOld;
    protected IAtom molecule;
    protected AtomSourceRandomMolecule moleculeSource;
    private AtomIteratorSinglet affectedAtomIterator;
    protected AtomGroupAction moveAction;
    protected IVector oldGeo, newGeo, temp, transVect;
    protected AtomPositionGeometricCenter centerer;
    protected MeterPotentialEnergy energyMeter;
    protected final IRandom random;
    
    public MCMoveCombinedCbmcTranslation(PotentialMaster pm, MCMoveCBMC mv,
            IRandom nRandom){
        super(pm);
        setCbmcMove(mv);
        this.random = nRandom;
        
        moleculeSource = new AtomSourceRandomMolecule();
        moleculeSource.setPhase(pm.getSimulation().getPhases()[0]);
        ((AtomSourceRandomMolecule)moleculeSource).setRandom(random);
        energyMeter = new MeterPotentialEnergy(pm);
        energyMeter.setPhase(pm.getSimulation().getPhases()[0]);
        
        affectedAtomIterator = new AtomIteratorSinglet();
        AtomActionTranslateBy translator = new AtomActionTranslateBy(pm.getSpace());
        transVect = translator.getTranslationVector();
        moveAction = new AtomGroupAction(translator);
        centerer = new AtomPositionGeometricCenter(pm.getSpace());
        
        oldGeo = pm.getSpace().makeVector();
        newGeo = pm.getSpace().makeVector();
        temp = pm.getSpace().makeVector();
    }
    public MCMoveCombinedCbmcTranslation(PotentialMaster pm, MCMoveCBMC mv, 
            IRandom nRandom, Phase ph){
        this(pm, mv, nRandom);
        setPhase(ph);
    }
    
    public AtomIterator affectedAtoms() { return affectedAtomIterator; }
    
    public double energyChange()  {
        return (uNew - uOld) + cbmcMove.energyChange();
    }

    public void acceptNotify() {
        // Nothing needs to be done!
    }

    public boolean doTrial() {
        molecule = moleculeSource.getAtom();
        affectedAtomIterator.setAtom(molecule);
        oldGeo.E(centerer.position(molecule));
        transVect.E(oldGeo);
        
        cbmcMove.doTrial(molecule);
        transVect.ME(centerer.position(molecule));

        uOld = energyMeter.getDataAsScalar();
        moveAction.actionPerformed(molecule);
        uNew = energyMeter.getDataAsScalar();
       
        return false;
    }

    public double getA() {
        return 1.0 * cbmcMove.getA();
    }

    public double getB() {
        return -(uNew - uOld) + cbmcMove.getB();
    }

    public void rejectNotify() {
        transVect.TE(-1.0);
        moveAction.actionPerformed(molecule);
        cbmcMove.rejectNotify();
    }

    public void setCbmcMove(MCMoveCBMC cbmcMove) {
        this.cbmcMove = cbmcMove;
    }

}
