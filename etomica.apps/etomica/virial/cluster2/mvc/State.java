package etomica.virial.cluster2.mvc;

import java.util.Set;

public interface State {

  public void clear();

  public Set<String> getKeys();

  public Object getProperty(String key);

  public void setProperty(String key, Object value);
}