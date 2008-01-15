package etomica.virial;

import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;
import etomica.species.Species;

public class SpeciesAlkane extends Species {

    public SpeciesAlkane(ISimulation sim, int numCarbons) {
        super();
        setMoleculeFactory(new AtomFactoryAlkane(this, sim.getSpace(), new ConformationLinear(sim)));
        AtomFactoryMono[] childFactories = new AtomFactoryMono[2];
        AtomTypeSphere atomTypeCH3 = new AtomTypeSphere(new ElementSimple("CH3", 15));
        atomTypeCH3.setParentType((AtomTypeMolecule)factory.getType());
        childFactories[0] = sim.isDynamic() ?
                   new AtomFactoryMonoDynamic(sim.getSpace(), atomTypeCH3) :
                   new AtomFactoryMono(sim.getSpace(), atomTypeCH3);
        AtomTypeSphere atomTypeCH2 = new AtomTypeSphere(new ElementSimple("CH2", 14));
        atomTypeCH2.setParentType((AtomTypeMolecule)factory.getType());
        childFactories[1] = sim.isDynamic() ?
                   new AtomFactoryMonoDynamic(sim.getSpace(), atomTypeCH2) :
                   new AtomFactoryMono(sim.getSpace(), atomTypeCH2);
        ((AtomFactoryAlkane)factory).setChildFactory(childFactories);
        ((AtomFactoryAlkane)factory).setTotalChildren(numCarbons);
        if (numCarbons > 1) {
            ((AtomFactoryAlkane)factory).setChildCount(new int[]{2,numCarbons-2});
        }
        else {
            ((AtomFactoryAlkane)factory).setChildCount(new int[]{1,0});
        }
    }
}
