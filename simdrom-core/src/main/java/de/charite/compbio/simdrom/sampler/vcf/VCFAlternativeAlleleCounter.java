/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import java.util.HashSet;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.filter.IFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFAlternativeAlleleCounter {

	private CloseableIterator<VariantContext> iterator;
	private ImmutableSet<IFilter> filters;
	private Optional<String> sample;
	private int counts = -1;

	public VCFAlternativeAlleleCounter(CloseableIterator<VariantContext> iterator, ImmutableSet<IFilter> filters,
			Optional<String> sample) {
		this.iterator = iterator;
		this.filters = filters;
		this.sample = sample;
	}

	/**
	 * @return The number of possible alternative alleles using all positions. If a sample is set only the genotype of
	 *         that sample is used.
	 */
	public int getCounts() {
		if (counts < 0)
			count();
		return counts;
	}

	private void count() {
		counts = 0;
		while (iterator.hasNext()) {
			Optional<VariantContext> optional_vc = Optional.of(iterator.next());
			for (IFilter iFilter : filters) {
				optional_vc = iFilter.filter(optional_vc);
			}
			if (optional_vc.isPresent()) {
				VariantContext vc = optional_vc.get();
				if (sample.isPresent()) {
					Set<Allele> alleleSet = new HashSet<>();
					alleleSet.addAll(vc.getGenotype(sample.get()).getAlleles());
					for (Allele allele : alleleSet) {
						if (allele.isNonReference())
							counts += 1;
					}
				} else
					counts += vc.getAlternateAlleles().size();
			}

		}
		// FIXME Really close it? test this!
		iterator.close();
	}
}
