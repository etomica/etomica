package etomica.units;

public abstract class BaseUnitPseudo3D extends BaseUnit implements BaseUnit.D3 {
    
    public static final double FALSE_DEPTH = 5.0;  //Angstroms

    private BaseUnit baseUnit3D;
    
    public BaseUnitPseudo3D(BaseUnit unit3D) {
        baseUnit3D = unit3D;
        name = unit3D.toString();
        symbol = unit3D.symbol();
        prefixAllowed = unit3D.prefixAllowed();
    }
    
    public Dimension dimension() {return baseUnit3D.dimension();}
    
    public double toPixels(double x) {return baseUnit3D.toPixels(x);}
    
    public static class Pressure extends BaseUnitPseudo3D {
        public Pressure(BaseUnit.Pressure baseUnit) {
            super(baseUnit);
            to = baseUnit.to * FALSE_DEPTH;
            from = 1.0/to;
        }
    }//end of Pressure
    
    public static class Volume extends BaseUnitPseudo3D {
        public Volume(BaseUnit.Volume baseUnit) {
            super(baseUnit);
            to = baseUnit.to / FALSE_DEPTH;
            from = 1.0/to;
        }
    }//end of Volume
    
}//end of BaseUnitPseudo3D