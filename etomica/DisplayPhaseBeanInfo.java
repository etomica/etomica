package etomica;
import etomica.beans.*;
import etomica.units.Unit;
import java.awt.*;
import java.beans.*;

public class DisplayPhaseBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Phase.class, PhaseEditor.class);
        //PropertyEditorManager.registerEditor(Unit.class, DeviceUnitEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("name",       DisplayPhase.class),
                new PropertyDescriptor("drawOverflow", DisplayPhase.class),
                new PropertyDescriptor("drawBoundary",      DisplayPhase.class),
                new PropertyDescriptor("drawMeters",      DisplayPhase.class),
                new PropertyDescriptor("drawSpace",       DisplayPhase.class),
                new PropertyDescriptor("imageShells",  DisplayPhase.class),
                new PropertyDescriptor("movable",       DisplayPhase.class),
                new PropertyDescriptor("resizable",       DisplayPhase.class),
                new PropertyDescriptor("scale",       DisplayPhase.class),
                new PropertyDescriptor("phase",       DisplayPhase.class),
                new PropertyDescriptor("updateInterval", DisplayPhase.class)};
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
    
    public Image getIcon(int iconType){
        String name = "";
        if(iconType == BeanInfo.ICON_COLOR_16x16){
            name="COLOR_16x16";
        }
        else if (iconType == BeanInfo.ICON_COLOR_32x32){
            name="COLOR_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_32x32){
            name="MONO_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_16x16){
            name="MONO_16x16";
        }
        else return null;
        return loadImage("Molecule_" + name + ".gif");
    }
}

