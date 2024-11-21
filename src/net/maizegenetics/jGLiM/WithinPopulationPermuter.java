package net.maizegenetics.jGLiM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

import net.maizegenetics.jGLiM.dm.FactorModelEffect;

public class WithinPopulationPermuter {
	final double[] data;
	final int npops;
	ArrayList<int[]> popIndices = new ArrayList<int[]>();
	static Random randomizer = new Random();
	
	public static void main(String[] args) {
		double[] od = new double[25];
		int[] pops = new int[25];
		int count = 0;
		for (int i = 0; i < 6; i++) {
			od[count] = 1.0 + ((double) i)/100.0;
			pops[count++] = 0;
		}
		for (int i = 0; i < 7; i++) {
			od[count] = 2.0 + ((double) i)/100.0;
			pops[count++] = 1;
		}
		for (int i = 0; i < 6; i++) {
			od[count] = 3.0 + ((double) i)/100.0;
			pops[count++] = 2;
		}
		for (int i = 6; i < 12; i++) {
			od[count] = 1.0 + ((double) i)/100.0;
			pops[count++] = 0;
		}
		

		System.out.print("Original:");
		for (int i = 0; i < 25; i++) System.out.format(" %1.2f", od[i]);
		System.out.println();
		
		FactorModelEffect fme = new FactorModelEffect(pops, true);
		WithinPopulationPermuter wpp = new WithinPopulationPermuter(od, fme);
		double[] perm = wpp.getPermutedData();
		
		System.out.print("Permuted:");
		for (int i = 0; i < 25; i++) System.out.format(" %1.2f", perm[i]);
		System.out.println();
		
		perm = wpp.getPermutedData();
		
		System.out.print("Permuted:");
		for (int i = 0; i < 25; i++) System.out.format(" %1.2f", perm[i]);
		System.out.println();

	}
	
	public WithinPopulationPermuter(double[] originalData, FactorModelEffect popEffect) {
		data = originalData;
		npops = popEffect.getNumberOfLevels();
		int[] levels = popEffect.getLevels();
		int[] levelCounts = popEffect.getLevelCounts();

		for(int p = 0; p < npops; p++) {
			popIndices.add(new int[levelCounts[p]]);
		}

		int[] count = new int[npops];

		int n = levels.length;
		for (int i = 0; i < n; i++) {
			int pop = levels[i];
			popIndices.get(pop)[count[pop]++] = i;
		}

	}
	
	public double[] getPermutedData() {
		int ndata = data.length;
		double[] permutedData = Arrays.copyOf(data, ndata);
		
		for (int p = 0; p < npops; p++) {
			int[] ndx = popIndices.get(p);
			
			int n = ndx.length;
			for (int i = n - 1; i >= 1; i--) {
				int j = randomizer.nextInt(i + 1);
				double temp = permutedData[ndx[j]];
				permutedData[ndx[j]] = permutedData[ndx[i]];
				permutedData[ndx[i]] = temp;
			}
		}
		
		return permutedData;
	}
	
	public void permuteData(double[] data) {
		int ndata = data.length;
		
		for (int p = 0; p < npops; p++) {
			int[] ndx = popIndices.get(p);
			
			int n = ndx.length;
			for (int i = n - 1; i >= 1; i--) {
				int j = randomizer.nextInt(i + 1);
				double temp = data[ndx[j]];
				data[ndx[j]] = data[ndx[i]];
				data[ndx[i]] = temp;
			}
		}
		
	}
}