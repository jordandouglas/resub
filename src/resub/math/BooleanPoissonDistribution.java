package resub.math;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;


@Description("Distribution of a boolean vector, where the number of trues comes from a Poisson distribution")
public class BooleanPoissonDistribution extends Distribution implements Loggable {
    final public Input<RealParameter> pInput = new Input<>("lambda", "mean number of 'true' values", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> xInput = new Input<>("x", "the boolean parameter", Input.Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        
    	
    	if (pInput.get().getValue() < 0) {
    		throw new IllegalArgumentException("lambda must be greater than 0");
    	}
    	
    	
    }
    
    
    @Override
    public double calculateLogP() {
        
    	double lambda = pInput.get().getValue();
    	
    	//org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(lambda);
    	
    	if (lambda < 0) {
    		logP = Double.NEGATIVE_INFINITY;
    		return logP;
    	}
    	
    	int count = getCount();
    	logP = count*Math.log(lambda) - lambda - logFactorial(count);
    	return logP;
    	
    }
    
    private double logFactorial(int n) {
    	double result = 0;
    	for (int i = 2; i <=n; i ++) {
    		result += Math.log(i);
    	}
    	return result;
    }

   
   
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(((BEASTObject)xInput.get()).getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(pInput.get().getID());
		return conds;
	}

	
	public int getCount() {
		int count = 0;
    	for (int i = 0; i <  xInput.get().getDimension(); i++) {
    		boolean val = xInput.get().getValue(i);
    		if (val) count++;
    	}
    	return count;
	}
	
	@Override
	public void init(PrintStream out) {
		out.print(xInput.get().getID() + ".prior\t" + xInput.get().getID() + ".count\t");
	}
	
	public void log(long sample, PrintStream out) {
		out.print(logP + "\t" + getCount() + "\t");
	}
	
	

	@Override
	public void sample(State state, Random random) {
		
		if (sampledFlag) return;
        sampledFlag = true;
        
        // Sample conditions
        sampleConditions(state, random);
        
        // Sample number of trues
        
        try {
        	org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(pInput.get().getValue());
			int k = dist.inverseCumulativeProbability(random.nextFloat());
			int n = xInput.get().getDimension();
			k = Math.min(k, n);
			
			// Randomly assign k trues
			List<Integer> indices = new ArrayList<>();
			for (int i = 0; i < k; i ++) {
				int index = random.nextInt(n);
				while (indices.contains(index)) {
					index = random.nextInt(n);
				}
				indices.add(index);
			}
			
			// Set values
			for (int i = 0; i < n; i ++) {
				if (indices.contains(i)) {
					xInput.get().setValue(i, true);
				}else {
					xInput.get().setValue(i, false);
				}
			
			}
			
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

} 










