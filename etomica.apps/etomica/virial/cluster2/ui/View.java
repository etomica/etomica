package etomica.virial.cluster2.ui;

import java.util.Map;

public interface View {

  public void activate();

  public void configure();
  
  public Map<String, Object> getData();
}
