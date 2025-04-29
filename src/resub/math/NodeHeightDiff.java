package resub.math;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.CalculationNode;

@Description("Average node height difference")
public class NodeHeightDiff extends CalculationNode implements Loggable, Function {
	
	
	public final Input<List<MRCAPrior>> mrcaPriorInput = new Input<>("clade", "a clade", new ArrayList<>());

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		
		
		double meanDiff = 0;
		int nclades = mrcaPriorInput.get().size();
		for (int i = 0; i < nclades; i ++) {
			for (int j = 0; j < nclades; j ++) {
				if (i == j) continue;
				
				double h1 = mrcaPriorInput.get().get(i).getCommonAncestor().getHeight();
				double h2 = mrcaPriorInput.get().get(j).getCommonAncestor().getHeight();
				double diff = Math.abs(h1 - h2);
				meanDiff += diff;
				
			}
		}
		
		meanDiff = meanDiff / (nclades*(nclades-1));
		
		
		return meanDiff;
	}

	@Override
	public void init(PrintStream out) {
		out.print(this.getID() + "\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(this.getArrayValue() + "\t");
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
