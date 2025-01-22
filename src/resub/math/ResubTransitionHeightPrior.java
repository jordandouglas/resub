package resub.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.evolution.tree.Node;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;



@Description("Places a Laplace prior with a mean of transition height")
public class ResubTransitionHeightPrior extends Distribution {
	
	final public Input<RealParameter> m_x = new Input<>("x", "point at which the density is calculated", Validate.REQUIRED);
	final public Input<RealParameter> rateInput = new Input<>("rate", "the rate of the Laplace", Validate.REQUIRED);
	final public Input<MRCAPrior> cladeInput = new Input<>("clade", "the clade with the node height", Validate.REQUIRED);
	
	
	@Override
    public void initAndValidate() {
		
		
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
    	
    	logP = 0;
    	
    	// Get common ancestor
    	Node mrca = cladeInput.get().getCommonAncestor();
    	
    	// Get height of this clade
    	double height = mrca.getHeight();
    	
    	
    	// Laplace dist on the distance between x and that height
    	double x = m_x.get().getValue();
    	double xadj = Math.abs(x - height);
    	double rate = rateInput.get().getValue();
    	
    	if (x <= 0 || x >= cladeInput.get().treeInput.get().getRoot().getHeight()) {
    		logP = Double.NEGATIVE_INFINITY;
    	}else {
    		logP = Math.log(rate) - Math.log(2) - rate*xadj;
    	}
    	
    	
    	
    	//if (x < height) {
    		//Log.warning(x + " / " + height + " -> " + xadj + " " + logP);
    	//}
    	return logP;
    	
    }
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(m_x.get().getID());
		return args;
	}
	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(rateInput.get().getID());
		conds.add(cladeInput.get().treeInput.get().getID());
		return conds;
	}
	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	   

}
