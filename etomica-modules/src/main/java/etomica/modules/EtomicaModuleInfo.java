package etomica.modules;

public final class EtomicaModuleInfo {
    public final String moduleName;
    public final String description;
    public final Class<?> moduleClass;
    public final String args;

    public EtomicaModuleInfo(String moduleName, String description, Class<?> moduleClass, String args) {
        this.moduleName = moduleName;
        this.description = description;
        this.moduleClass = moduleClass;
        this.args = args;
    }

    public static EtomicaModuleInfo mod(String moduleName, String description, Class<?> moduleClass) {
        return new EtomicaModuleInfo(moduleName, description, moduleClass, "");
    }

    public static EtomicaModuleInfo mod(String moduleName, String description, Class<?> moduleClass, String args) {
        return new EtomicaModuleInfo(moduleName, description, moduleClass, args);
    }

    public String toString() {
        return this.moduleName;
    }

    public static final EtomicaModuleInfo[] ETOMICA_MODULES = {
            mod("Adsorption", "No description yet", etomica.modules.adsorption.AdsorptionFastererGraphic.class),
            mod("Catalysis", "No description yet", etomica.modules.catalysis.CatalysisGraphicFasterer.class),
            mod("Chain Equilibrium", "No description yet", etomica.modules.chainequilibrium.ChainEquilibriumFastererGraphic.class),
            mod("Colloid", "No description yet", etomica.modules.colloid.ColloidGraphicFasterer.class),
            mod("Crystal Viewer", "Permits viewing of static 3-dimensional lattices and the planes they define. Various elementary lattices can be selected, and any plane or surface can be viewed by specifying it via its Miller indices. Image can be rotated to permit viewing from any angle.", etomica.modules.crystalviewer.CrystalViewer.class),
            mod("Discontinuous MD 2D", "No description yet", etomica.modules.swmd.SwmdGraphicFasterer.class, "2"),
            mod("Discontinuous MD 3D", "No description yet", etomica.modules.swmd.SwmdGraphicFasterer.class, "3"),
            mod("DCVGCMD", "No description yet", etomica.modules.dcvgcmd.DCVGCMDGraphicFasterer.class),
            mod("Droplet", "No description yet", etomica.modules.droplet.DropletGraphicFasterer.class),
            mod("Droplet (atomic)", "No description yet", etomica.modules.droplet.DropletAtomicGraphicFasterer.class),
            mod("Ensembles", "No description yet", etomica.modules.ensembles.LJMCGraphicFasterer.class),
            mod("Glass", "No description yet", etomica.modules.glass.GlassGraphicFasterer.class),
            mod("Insertion", "No description yet", etomica.modules.insertion.InsertionGraphicFasterer.class),
            mod("Interfacial", "No description yet", etomica.modules.interfacial.InterfacialSWGraphicFasterer.class),
            mod("LJMD", "No description yet", etomica.modules.ljmd.LjmdGraphicFasterer.class),
            mod("Material Fracture", "No description yet", etomica.modules.materialfracture.MaterialFractureGraphicFasterer.class),
            mod("Chemical Potential", "No description yet", etomica.modules.mu.MuGraphicFasterer.class),
            mod("Multiharmonic", "No description yet", etomica.modules.multiharmonic.MultiharmonicGraphicMCFasterer.class),
            mod("Nucleation", "No description yet", etomica.modules.nucleation.NucleationGraphic.class),
            mod("Osmosis", "No description yet", etomica.modules.osmosis.OsmosisFasterer.class),
            mod("Piston-Cylinder 2D", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a chamber with a movable wall under external pressure.", etomica.modules.pistoncylinder.PistonCylinderGraphicFasterer.class, "-dim 2"),
            mod("Piston-Cylinder 3D", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a 3D chamber with a movable wall under external pressure.", etomica.modules.pistoncylinder.PistonCylinderGraphicFasterer.class, "-dim 3"),
            mod("Piston-Cylinder-SWMD", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a chamber with a movable wall under external pressure.  This version of the simulation has extra features used for the SWMD module", etomica.modules.pistoncylinder.PistonCylinderGraphicFasterer.class, "-dim 2 -configButton -densityInput -nMoleculeSlider -rdf -fastButton"),
            mod("Reaction Equilibrium", "Simple reaction equlibrium involving two atomic species and the three dimeric molecules they can form. Atoms move about via 2-D molecular dynamics, and can \"react\" to form dimers. Dynamic equilibria is demonstrated through constant recombining of atoms, and equilibria can be quantified and analyzed with thermodynamic reaction equilibria models.", etomica.modules.reactionequilibrium.ReactionEquilibriumGraphicFasterer.class),
            mod("Rheology", "No description yet", etomica.modules.rheology.RheologyGraphic.class),
            mod("Reverse Osmosis", "No description yet", etomica.modules.rosmosis.ReverseOsmosisGraphicFasterer.class),
            mod("Reverse Osmosis (water)", "No description yet", etomica.modules.rosmosis.ReverseOsmosisWaterGraphicFasterer.class),
            mod("Self-Assembled Monolayer", "No description yet", etomica.modules.sam.SamGraphicFasterer.class),
            mod("Statistics", "No description yet", etomica.modules.statistics.StatisticsMCGraphicFasterer.class, "-moduleNum 1"),
            mod("Statistics Part 2", "No description yet", etomica.modules.statistics.StatisticsMCGraphicFasterer.class, "-moduleNum 2"),
            mod("VLE", "No description yet", etomica.modules.vle.VLEFasterer.class),
            mod("B2 fit (for VLE)", "No description yet", etomica.modules.vle.B2Fit.class),
    };
}

