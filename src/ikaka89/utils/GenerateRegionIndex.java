package ikaka89.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * 
 * @author akakumanu
 *
 */
public class GenerateRegionIndex {
	
	public int binwidth = 200;
	public int binstep = 200;
	
	public GenomeConfig gcon;
	
	public void setBinWidth(int bw) {binwidth = bw;}
	public void setBinStep(int bs){binstep=bs;}
	
	public GenerateRegionIndex(GenomeConfig gc) {
		gcon=gc;
	}
	
	
	public void executeByChrom() {
		Map<String,Integer> chrom_lengths = gcon.getGenome().getChromLengthMap();
		int index = 0;
		for(String chr : chrom_lengths.keySet()) {
			int length = chrom_lengths.get(chr);
			for(int midpt=binwidth/2; midpt<length; midpt=midpt+binstep) {
				Region reg = new Region(gcon.getGenome(), chr, midpt-binwidth/2, midpt+binwidth/2);
				index++;
				StringBuilder  sb = new StringBuilder();
				sb.append(reg.toString());sb.append("\t");sb.append(index);
				System.out.println(sb.toString());
			}
		}
	}
	
	public static void main(String[] args) {
		GenomeConfig gc = new GenomeConfig(args);
		
		GenerateRegionIndex runner = new GenerateRegionIndex(gc);
		runner.setBinStep(Args.parseInteger(args, "binstep", 200));
		runner.setBinWidth(Args.parseInteger(args, "binwidth", 200));
		
		runner.executeByChrom();
		
	}
	
	
}
