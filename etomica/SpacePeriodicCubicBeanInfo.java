package simulate;
import java.awt.*;
import java.beans.*;



public class SpacePeriodicCubicBeanInfo extends SimpleBeanInfo {
    public Image getIcon(int iconType){
       // String name = "";
        //if(iconType == BeanInfo.ICON_COLOR_16X16){
        //    name="COLOR_16X16";
        //}
        //else if (iconType == BeanInfo.ICON_COLOR_32X32){
        //    name="COLOR_32X32";
       // }
       // else if (iconType == BeanInfo.ICON_MONO_32X32){
       //     name="MONO_32X32";
       // }
        //else if (iconType == BeanInfo.ICON_MONO_16X16){
         //   name="MONO_16X16";
       // }
       // else return null;
        return loadImage("sim.gif");
    }
}
