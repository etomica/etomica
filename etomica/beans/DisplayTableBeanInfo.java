package etomica;
import etomica.beans.*;
import etomica.units.Unit;
import java.awt.*;
import java.beans.*;

public class DisplayTableBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Meter[].class, MeterArrayEditor.class);
        PropertyEditorManager.registerEditor(MeterFunction.class, MeterFunctionEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("name",       DisplayTable.class),
                new PropertyDescriptor("meters", DisplayTable.class),
                new PropertyDescriptor("meterFunction", DisplayTable.class),
                new PropertyDescriptor("resetVisible",      DisplayTable.class),
                new PropertyDescriptor("showingKE",      DisplayTable.class),
                new PropertyDescriptor("showingPE",       DisplayTable.class),
                new PropertyDescriptor("label",  DisplayTable.class),
                new PropertyDescriptor("updateInterval", DisplayTable.class)};
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

