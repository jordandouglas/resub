package resub.math;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.CalculationNode;

@Description("Difference in node height")
public class NodeHeightDiff extends CalculationNode implements Loggable, Function {
	
	
	public final Input<MRCAPrior> mrcaPrior1Input = new Input<>("clade1", "the clade", Validate.REQUIRED);
	public final Input<MRCAPrior> mrcaPrior2Input = new Input<>("clade2", "the clade", Validate.REQUIRED);

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
		double h1 = mrcaPrior1Input.get().getCommonAncestor().getHeight();
		double h2 = mrcaPrior2Input.get().getCommonAncestor().getHeight();
		return h1 - h2;
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
