package etomica;

public interface Action {
    
    public void actionPerformed();
    
    public interface Undoable {
        public void attempt();
        public void undo();
    }
    
    public String getLabel();
    public void setLabel(String text);
}