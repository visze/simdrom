/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class VCFSampler implements Iterator<VariantContext> {

	private double probability;
	private List<Integer> selectAlleles;
	private int position = -1;
	private String afIdentifier;
	private String sample = null;
	private VCFFileReader parser;
	private CloseableIterator<VariantContext> iterator;
	private Random random;
	private String filePath;

	public VCFSampler(String filePath) {
		this.filePath = filePath;
		this.parser = new VCFFileReader(new File(filePath), false);
	}

	public CloseableIterator<VariantContext> getIterator() {
		if (this.iterator == null)
			this.iterator = this.parser.iterator();
		return iterator;
	}

	@Override
	public boolean hasNext() {
		// FIXME has next can be true, but next can gi9ve back null!
		return getIterator().hasNext();
	}

	@Override
	public VariantContext next() {
		return getNextVariant();
	}

	private VariantContext getNextVariant() {
		VariantContext output = null;
		while (getIterator().hasNext() && output == null) {
			VariantContext candidate = getIterator().next();
			Map<Integer, Boolean> alleles = useAlleles(candidate);
			if (!alleles.isEmpty()) {
				output = createVariantContextWithGenotype(candidate, alleles);
				break;
			}
		}
		return output;
	}

	private VariantContext createVariantContextWithGenotype(VariantContext candidate, Map<Integer, Boolean> alleles) {
		if (useSample()) {
			if (candidate.hasGenotype(getSample()))
				return new VariantContextBuilder(candidate).alleles(getAlleles(candidate, alleles.keySet()))
						.genotypes(candidate.getGenotypes(getSample())).make();
			else
				return null;
		} else {
			return new VariantContextBuilder(candidate).genotypes(createGenotype(candidate.getAlleles(),alleles)).make();
		}
	}

	private Collection<Allele> getAlleles(VariantContext candidate, Set<Integer> posOfAllele) {
		Collection<Allele> alleles = new ArrayList<Allele>();
		alleles.add(candidate.getReference());
		for (int i : posOfAllele) {
			alleles.add(candidate.getAlternateAlleles().get(i));
		}
		return alleles;
	}

	private Genotype createGenotype(List<Allele> alleles, Map<Integer, Boolean> use) {
		List<Allele> filteredAlleles = new ArrayList<Allele>();

		int allele = use.keySet().iterator().next() +1;
		// more the one alternative allele, do not use ref!
		if (use.size() > 1) {
			for (int i : use.keySet()) {
				filteredAlleles.add(alleles.get(i+1));
			}
		} else if (use.get(allele-1)){
			
			filteredAlleles.add(alleles.get(allele));
			filteredAlleles.add(alleles.get(allele));
		} else {
			filteredAlleles.add(alleles.get(0));
			filteredAlleles.add(alleles.get(allele));
			
		}

		return GenotypeBuilder.create("Sampled", filteredAlleles);
	}

	private Map<Integer, Boolean> useAlleles(VariantContext candidate) {
		Map<Integer, Boolean> candidates = new HashMap<Integer, Boolean>();

		if (useAF()) {// AF flag
			Object af = candidate.getCommonInfo().getAttribute(getAFIdentifier());
			if (af instanceof ArrayList<?>) {
				if (((ArrayList<?>) af).get(0) instanceof String) {
					int i = 0;
					for (Object o : (ArrayList<?>) af) {
						addCandidateByHardyWeinberg(candidates, i, Double.parseDouble((String) o));
						i++;
					}
				}
			} else {
				addCandidateByHardyWeinberg(candidates, 0,
						candidate.getCommonInfo().getAttributeAsDouble(getAFIdentifier(), 0.0));
			}
		} else if (useCounts()) { // probability > 1
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				this.position++;
				if (selectAlleles.contains(position + i))
					candidates.put(i, nextDouble() <= 0.5);
			}

		} else { // probability
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				addCandidateByHardyWeinberg(candidates, 0, getProbability());
			}
		}
		return candidates;

	}

	private void addCandidateByHardyWeinberg(Map<Integer, Boolean> candidates, int i, double af) {
		double[] homHetAF = getHardyWeinbergPrincipleHomhet(af);
		double random = nextDouble();
		if (random <= homHetAF[0])
			candidates.put(i, true);
		else if (random <= homHetAF[1])
			candidates.put(i, false);
	}

	private double[] getHardyWeinbergPrincipleHomhet(double af) {
		return new double[]{Math.pow(1.0-Math.sqrt(1.0-af),2),af};
	}

	private boolean useCounts() {
		return probability > 1.0;
	}

	private double nextDouble() {
		if (random == null)
			random = new Random();
		return random.nextDouble();
	}

	public void setProbability(double probability) {
		this.probability = probability;
		if (useCounts()) {
			VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(filePath);
			setCounts(counter.getCounts());
		}
	}

	public void setSample(String sample) {
		this.sample = sample;
	}

	public String getSample() {
		return sample;
	}

	public double getProbability() {
		return probability;
	}

	public VCFHeader getFileHeader() {
		// List<VCFHeaderLine> headerLines = new ArrayList<VCFHeaderLine>();
		// headerLines.addAll();

		Set<VCFHeaderLine> set = new LinkedHashSet<VCFHeaderLine>();
		set.addAll(parser.getFileHeader().getMetaDataInInputOrder());
		set.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		return new VCFHeader(set, getSampleNames());
	}

	public ImmutableSet<String> getSampleNames() {
		if (getSample() == null)
			return ImmutableSet.<String> builder().add("Sampled").build();
		return ImmutableSet.<String> builder().add(getSample()).build();
	}

	private void setCounts(int counts) {
		List<Integer> randomAlleles = new ArrayList<Integer>(counts);
		for (int i = 0; i < counts; i++) {
			randomAlleles.add(i);
		}
		Collections.shuffle(randomAlleles);
		this.selectAlleles = new ArrayList<Integer>((int) Math.floor(getProbability()));
		for (int i = 0; i < (int) Math.floor(getProbability()); i++) {
			this.selectAlleles.add(randomAlleles.get(i));
		}
		Collections.sort(this.selectAlleles);
	}

	public void setAFIdentifier(String afIdentifier) {
		this.afIdentifier = afIdentifier;
	}

	public String getAFIdentifier() {
		return afIdentifier;
	}

	private boolean useAF() {
		return getAFIdentifier() != null;
	}

	private boolean useSample() {
		return getSample() != null;
	}

	public void close() {
		parser.close();
	}

}
