package gld.algo.tlc;

public class TrafficState {
	
    public int[] qLengths;  

    /**
     * Creates a new <code>GridState</code> instance with the default
     * start row and column. */
    public TrafficState(int[] x) 
    {
	    qLengths = new int[x.length];
	    for (int i = 0; i < x.length; i++) 
	    {
			qLengths[i] = x[i];
		}
    }
    
    public boolean equals(TrafficState state){
    	boolean isEqual = true;
    	
    	for (int i = 0; i < qLengths.length; i++) 
    	{
    		if(qLengths[i] != state.qLengths[i]) {isEqual=false;break;}			
		}
    	return isEqual;
        }   
}
