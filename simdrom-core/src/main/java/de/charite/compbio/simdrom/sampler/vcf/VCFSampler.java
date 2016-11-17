/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.filter.IFilter;
import de.charite.compbio.simdrom.sampler.DeNovoSampler;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeEOFException;
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

/**
 * 
 * The {@link VCFSampler} class uses a VCF file as basis and acts as an iterator over such a file. There are multiple
 * possible options what the {@link VCFSampler} does with the read {@link VariantContext}s.
 * <p>
 * <dl>
 * <dt>Filter</dt>
 * <dd>Variants can be filtered using multiple {@link IFilter}. Completely removed variants where not finally passed
 * into the iterator.</dd>
 * <dt>Sample selection</dt>
 * <dd>If a multi-VCF file is used there is the possibility to select one sample out of it. It can be done randomly by
 * running the VCFRandomSampleSelecter or you can predefine a sample name. The sample name must be present in the VCF
 * header.</dd>
 * <dt>Variant generation using Hardy-Weinberg principle</dt>
 * <dd>If allele frequencies are present in the pares VCF they can be used to generate a completely new artificial
 * {@link Genotype} by using the Hardy-Weinberg principle. Direct allele frequencies can be used of the allele count
 * divided by the number of total alleles.</dd> *
 * <dt>Interval list</dt>
 * <dd>If an {@link IntervalList} is present only variants of that lists will be pares by the {@link VCFFileReader} and
 * these variants are only used by the {@link VCFSampler}.</dd>
 * <dt>De-novo generation</dt>
 * <dd>Inventing new variations. This function is not implemented yet.</dd>
 * </dl>
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFSampler implements CloseableIterator<VariantContext> {

	/**
	 * The probability to select a variant. If set to 1 all variants will be selected.
	 */
	private double probability;
	/**
	 * 
	 */
	private List<Integer> selectAlleles;
	/**
	 * The maximum number of variants that should be selected by the sampler. If set to smaller than 1 this option is of
	 * and there will be no maximum limit used.
	 */
	private int variantsAmount;
	private int position = -1;
	/**
	 * The identifier of the allele frequency in the Info column.
	 */
	private String afIdentifier;
	/**
	 * The identifier of the alt allele counts in the Info column.
	 */
	private String acIdentifier;
	/**
	 * The identifier of the total number of alleles.
	 */
	private String anIdentifier;
	/**
	 * If set the genotypes of this samples are selected. if null no sample is selected.
	 */
	private String sample;
	/**
	 * The VCF file reader.
	 */
	private VCFFileReader reader;

	/**
	 * The variant iterator defined by the reader and maybe by intervals.
	 */
	private CloseableIterator<VariantContext> iterator;
	/**
	 * The next variantContext. If this is null there will be no further variant and {@link #hasNext()} is false.
	 */
	private VariantContext next;
	private Random random;
	private DeNovoSampler deNovoGenerator;
	private ImmutableSet<IFilter> filters;
	// intervals
	private IntervalList intervals;
	private int intervalPosition = 0;

	/**
	 * Default constructor initializing all fields. Can only be used by the {@link Builder}.
	 * 
	 * @param reader
	 *            {@link #probability}
	 * @param probability
	 *            {@link #probability}
	 * @param sample
	 *            {@link #sample}
	 * @param variantsAmount
	 *            {@link #variantsAmount}
	 * @param afIdentifier
	 *            {@link #afIdentifier}
	 * @param acIdentifier
	 *            {@link #acIdentifier}
	 * @param anIdentifier
	 *            {@link #anIdentifier}
	 * @param filters
	 *            {@link #filters}
	 * @param intervals
	 *            {@link #intervals}
	 * @param deNovoSampler
	 *            {@link #deNovoGenerator}
	 * @param selectAlleles
	 *            {@link #selectAlleles}
	 * @param random
	 *            {@link #random}
	 */
	private VCFSampler(VCFFileReader reader, double probability, String sample, int variantsAmount, String afIdentifier,
			String acIdentifier, String anIdentifier, ImmutableSet<IFilter> filters, IntervalList intervals,
			DeNovoSampler deNovoSampler, List<Integer> selectAlleles, Random random) {
		this.reader = reader;
		this.probability = probability;
		this.sample = sample;
		this.variantsAmount = variantsAmount;
		this.afIdentifier = afIdentifier;
		this.acIdentifier = acIdentifier;
		this.anIdentifier = anIdentifier;
		this.filters = filters;
		this.intervals = intervals;
		this.selectAlleles = selectAlleles;
		this.deNovoGenerator = deNovoSampler;
		this.random = random;
		this.next = getNextVariant();
	}

	public static final class Builder {

		private VCFFileReader reader;
		private double probability = 1.0;
		private List<Integer> selectAlleles;
		private String sample = null;
		private int variantsAmount = 0;
		private String afIdentifier = null;
		private String acIdentifier = null;
		private String anIdentifier = null;
		private ImmutableSet<IFilter> filters = ImmutableSet.<IFilter> builder().build();
		IntervalList intervals = new IntervalList(new SAMFileHeader());
		DeNovoSampler deNovoSampler;
		private long seed;
		private boolean useSeed = false;

		public Builder() {
		}

		public Builder file(String path) {
			this.reader = new VCFFileReader(new File(path));
			return this;
		}

		public Builder file(File file) {
			this.reader = new VCFFileReader(file);
			return this;
		}

		public Builder vcfReader(VCFFileReader reader) {
			this.reader = reader;
			return this;
		}

		public Builder probability(double probability) {
			this.probability = probability;
			return this;
		}

		public Builder sample(String sample) {
			this.sample = sample;
			return this;
		}

		public Builder afIdentifier(String afIdentifier) {
			this.afIdentifier = afIdentifier;
			return this;
		}

		public Builder acIdentifier(String acIdentifier) {
			this.acIdentifier = acIdentifier;
			return this;
		}

		public Builder anIdentifier(String anIdentifier) {
			this.anIdentifier = anIdentifier;
			return this;
		}

		public Builder variantsAmount(int variantsAmount) {
			this.variantsAmount = variantsAmount;
			return this;
		}

		public Builder intervals(IntervalList intervals) {
			this.intervals = intervals.uniqued().sorted();
			return this;
		}

		public Builder deNovoGenerator(DeNovoSampler deNovoSampler) {
			this.deNovoSampler = deNovoSampler;
			return this;
		}

		public Builder filters(ImmutableSet<IFilter> filters) {
			this.filters = filters;
			return this;
		}

		public Builder seed(long seed) {
			this.seed = seed;
			this.useSeed = true;
			return this;
		}

		public VCFSampler build() {

			// returns null if the reader is not set (not useful)
			if (reader == null)
				throw new RuntimeEOFException("No variants are set for the sampler: The reader cannot be null");

			// do not use it if one of them is not set!
			if (acIdentifier == null || anIdentifier == null) {
				acIdentifier = null;
				anIdentifier = null;
			}

			if (variantsAmount < 0)
				variantsAmount = 0;
			else if (variantsAmount > 0) {
				VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(reader.iterator(), filters);
				setCounts(counter.getCounts());
			}

			Random random;
			if (useSeed)
				random = new Random(seed);
			else
				random = new Random();

			return new VCFSampler(reader, probability, sample, variantsAmount, afIdentifier, acIdentifier, anIdentifier,
					filters, intervals, deNovoSampler, selectAlleles, random);
		}

		private void setCounts(int counts) {
			List<Integer> randomAlleles = new ArrayList<Integer>(counts);
			for (int i = 0; i < counts; i++) {
				randomAlleles.add(i);
			}
			Collections.shuffle(randomAlleles);
			this.selectAlleles = new ArrayList<Integer>(this.variantsAmount);
			for (int i = 0; i < this.variantsAmount; i++) {
				this.selectAlleles.add(randomAlleles.get(i));
			}
			Collections.sort(this.selectAlleles);
		}

	}

	public CloseableIterator<VariantContext> getIterator() {

		if (this.iterator == null) {
			if (useIntervals())
				this.iterator = getNextIntervalInterator();
			else
				this.iterator = this.reader.iterator();
		}
		return iterator;
	}

	private CloseableIterator<VariantContext> getNextIntervalInterator() {
		Interval interval = nextInterval();
		if (interval != null)
			return this.reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
		else
			return null;

	}

	private Interval nextInterval() {
		Interval output = null;
		if (intervalPosition < getIntervals().size()) {
			output = getIntervals().getIntervals().get(intervalPosition);
			intervalPosition++;
		}
		return output;
	}

	@Override
	public boolean hasNext() {
		return next != null;
	}

	/**
	 * Has next can be true, but next can give back null! But it is important to set the next iterator.
	 * 
	 * @return boolean
	 */
	private boolean checkForNext() {
		// FIXME
		if (useIntervals())
			while (getIterator() != null && !getIterator().hasNext()) {
				// close the actual
				this.iterator.close();
				// get a new one
				this.iterator = getNextIntervalInterator();
			}
		if (getIterator() == null)
			return false;
		return getIterator().hasNext();
	}

	private boolean useIntervals() {
		return !getIntervals().getIntervals().isEmpty();
	}

	@Override
	public VariantContext next() {
		VariantContext actual = next;
		next = getNextVariant();
		return actual;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	private VariantContext getNextVariant() {
		VariantContext output = null;
		while (checkForNext() && output == null) {
			// get next line
			VariantContext candidate = getIterator().next();

			// filter
			candidate = filter(candidate);
			if (candidate == null)
				continue;

			// get alleles by sampling method
			Map<Integer, Boolean> alleles = useAlleles(candidate);
			if (!alleles.isEmpty()) {
				output = createVariantContextWithGenotype(candidate, alleles);
				if (output == null)
					continue;
				break;
			}
		}
		return output;
	}

	private VariantContext filter(VariantContext candidate) {
		for (IFilter iFilter : getFilters()) {
			candidate = iFilter.filter(candidate);
		}
		return candidate;
	}

	private VariantContext createVariantContextWithGenotype(VariantContext candidate, Map<Integer, Boolean> alleles) {
		if (useSample()) {
			Genotype genotype = candidate.getGenotype(getSample());
			if (!genotype.isHomRef())
				return new VariantContextBuilder(candidate).genotypes(candidate.getGenotypes(getSample())).make();
			else
				return null;
		} else {
			return new VariantContextBuilder(candidate).genotypes(createGenotype(candidate.getAlleles(), alleles))
					.make();
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

		int allele = use.keySet().iterator().next() + 1;
		// more the one alternative allele, do not use ref!
		if (use.size() > 1) {
			for (int i : use.keySet()) {
				filteredAlleles.add(alleles.get(i + 1));
			}
		} else if (use.get(allele - 1)) {

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
		} else if (useAC()) {
			Object ac = candidate.getCommonInfo().getAttribute(getACIdentifier());
			int an = candidate.getCommonInfo().getAttributeAsInt(getANIdentifier(), 0);
			if (ac instanceof ArrayList<?>) {
				if (((ArrayList<?>) ac).get(0) instanceof String) {
					int i = 0;
					for (Object o : (ArrayList<?>) ac) {
						addCandidateByHardyWeinberg(candidates, i, Double.parseDouble((String) o) / (double) an);
						i++;
					}
				}
			} else {
				addCandidateByHardyWeinberg(candidates, 0,
						(double) candidate.getCommonInfo().getAttributeAsInt(getACIdentifier(), 0) / (double) an);
			}

		} else if (useCounts()) { // variantsAmount > 0
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

	private boolean useAC() {
		return getACIdentifier() != null && getANIdentifier() != null;
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
		return new double[] { Math.pow(1.0 - Math.sqrt(1.0 - af), 2), af };
	}

	private boolean useCounts() {
		return getVariantsAmount() > 0;
	}

	private double nextDouble() {
		return random.nextDouble();
	}

	public String getSample() {
		return sample;
	}

	public double getProbability() {
		return probability;
	}

	public VCFHeader getFileHeader() {

		Set<VCFHeaderLine> set = new LinkedHashSet<VCFHeaderLine>();
		set.addAll(reader.getFileHeader().getMetaDataInInputOrder());
		set.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		return new VCFHeader(set, getSampleNames());
	}

	public ImmutableSet<String> getSampleNames() {
		if (getSample() == null)
			return ImmutableSet.<String> builder().add("Sampled").build();
		return ImmutableSet.<String> builder().add(getSample()).build();
	}

	public int getVariantsAmount() {
		return variantsAmount;
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
		iterator.close();
		reader.close();
	}

	public ImmutableSet<IFilter> getFilters() {
		return filters;
	}

	public String getACIdentifier() {
		return acIdentifier;
	}

	public String getANIdentifier() {
		return anIdentifier;
	}

	public IntervalList getIntervals() {
		return intervals;
	}

}
